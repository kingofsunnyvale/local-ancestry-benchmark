"""Evaluation metrics for local ancestry inference benchmarking."""

import pandas as pd

import config
from src.parsers import (
    get_gnomix_msp_path,
    get_rfmix_msp_path,
    intervals_to_df,
    msp_to_intervals,
    parse_msp,
)


def _overlap_length(s1: int, e1: int, s2: int, e2: int) -> int:
    """Compute overlap length between two intervals."""
    return max(0, min(e1, e2) - max(s1, s2))


def per_base_accuracy(
    pred_df: pd.DataFrame,
    truth_df: pd.DataFrame,
) -> dict:
    """Compute per-base accuracy between predicted and ground-truth intervals.

    Both DataFrames must have columns: individual, haplotype, start, end, ancestry.

    Returns dict with overall accuracy and per-ancestry accuracy.
    """
    total_bases = 0
    correct_bases = 0
    per_anc_total = {}
    per_anc_correct = {}

    # Group both by (individual, haplotype) for efficient comparison
    pred_grouped = pred_df.groupby(["individual", "haplotype"])
    truth_grouped = truth_df.groupby(["individual", "haplotype"])

    for key, truth_group in truth_grouped:
        if key not in pred_grouped.groups:
            # No predictions for this haplotype — count all as incorrect
            for _, trow in truth_group.iterrows():
                span = trow["end"] - trow["start"]
                total_bases += span
                anc = trow["ancestry"]
                per_anc_total[anc] = per_anc_total.get(anc, 0) + span
            continue

        pred_group = pred_grouped.get_group(key)
        pred_list = list(pred_group.itertuples(index=False))

        for _, trow in truth_group.iterrows():
            t_start, t_end, t_anc = trow["start"], trow["end"], trow["ancestry"]
            span = t_end - t_start
            total_bases += span
            per_anc_total[t_anc] = per_anc_total.get(t_anc, 0) + span

            for prow in pred_list:
                overlap = _overlap_length(t_start, t_end, prow.start, prow.end)
                if overlap > 0 and prow.ancestry == t_anc:
                    correct_bases += overlap
                    per_anc_correct[t_anc] = per_anc_correct.get(t_anc, 0) + overlap

    overall_acc = correct_bases / total_bases if total_bases > 0 else 0.0

    per_anc_acc = {}
    for anc in sorted(per_anc_total.keys()):
        total = per_anc_total[anc]
        correct = per_anc_correct.get(anc, 0)
        per_anc_acc[anc] = correct / total if total > 0 else 0.0

    return {
        "overall": overall_acc,
        "per_ancestry": per_anc_acc,
        "total_bases": total_bases,
        "correct_bases": correct_bases,
    }


def concordance(
    pred_a_df: pd.DataFrame,
    pred_b_df: pd.DataFrame,
) -> dict:
    """Compute per-base concordance between two sets of predictions.

    Same algorithm as per_base_accuracy but comparing tool A vs tool B
    (neither is ground truth).
    """
    # Use pred_a as "truth" and pred_b as "prediction"
    return per_base_accuracy(pred_b_df, pred_a_df)


def global_ancestry_fractions(pred_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Compute global ancestry fractions per individual.

    Returns {individual: {ancestry: fraction}}.
    """
    results = {}
    for indiv, group in pred_df.groupby("individual"):
        total_len = 0
        anc_len = {}
        for _, row in group.iterrows():
            span = row["end"] - row["start"]
            total_len += span
            anc = row["ancestry"]
            anc_len[anc] = anc_len.get(anc, 0) + span

        fracs = (
            {anc: length / total_len for anc, length in anc_len.items()}
            if total_len > 0
            else {}
        )
        results[indiv] = fracs

    return results


def mean_global_fractions(pred_df: pd.DataFrame) -> dict[str, float]:
    """Compute mean global ancestry fractions across all individuals."""
    per_indiv = global_ancestry_fractions(pred_df)
    if not per_indiv:
        return {}

    all_ancs = set()
    for fracs in per_indiv.values():
        all_ancs.update(fracs.keys())

    mean_fracs = {}
    n = len(per_indiv)
    for anc in sorted(all_ancs):
        mean_fracs[anc] = sum(fracs.get(anc, 0.0) for fracs in per_indiv.values()) / n

    return mean_fracs


def evaluate_simulated(chrom: int) -> pd.DataFrame:
    """Evaluate both tools on all simulated conditions.

    Returns a DataFrame with columns:
    tool, scenario, generation, overall_accuracy, + per-ancestry accuracy columns.
    """
    rows = []

    for scenario in config.SCENARIOS:
        for gen in config.ADMIXTURE_GENERATIONS:
            gt_path = (
                config.SIMULATED_DIR / scenario / f"gen_{gen}" / "ground_truth.csv"
            )
            if not gt_path.exists():
                continue
            truth_df = pd.read_csv(gt_path)

            for tool_name in ["rfmix", "gnomix"]:
                if tool_name == "rfmix":
                    msp_path = get_rfmix_msp_path(scenario, gen, chrom)
                else:
                    msp_path = get_gnomix_msp_path(scenario, gen)

                if not msp_path.exists():
                    continue

                result = parse_msp(msp_path, tool_name)
                pred_intervals = msp_to_intervals(result)
                pred_df = intervals_to_df(pred_intervals)

                metrics = per_base_accuracy(pred_df, truth_df)

                row = {
                    "tool": tool_name,
                    "scenario": scenario,
                    "generation": gen,
                    "overall_accuracy": round(metrics["overall"], 4),
                }
                for anc, acc in metrics["per_ancestry"].items():
                    row[f"accuracy_{anc}"] = round(acc, 4)

                rows.append(row)
                print(
                    f"  {tool_name} {scenario}/gen_{gen}: "
                    f"overall={metrics['overall']:.2%}"
                )

    df = pd.DataFrame(rows)
    return df


def evaluate_concordance(chrom: int) -> pd.DataFrame:
    """Compute concordance between RFMix and Gnomix on real data.

    Returns DataFrame with per-scenario concordance.
    """
    rows = []

    for scenario in config.SCENARIOS:
        rfmix_path = get_rfmix_msp_path(scenario, None, chrom)
        gnomix_path = get_gnomix_msp_path(scenario, None)

        if not rfmix_path.exists() or not gnomix_path.exists():
            continue

        rfmix_result = parse_msp(rfmix_path, "rfmix")
        gnomix_result = parse_msp(gnomix_path, "gnomix")

        rfmix_intervals = msp_to_intervals(rfmix_result)
        gnomix_intervals = msp_to_intervals(gnomix_result)

        rfmix_df = intervals_to_df(rfmix_intervals)
        gnomix_df = intervals_to_df(gnomix_intervals)

        metrics = concordance(rfmix_df, gnomix_df)

        row = {
            "scenario": scenario,
            "concordance": round(metrics["overall"], 4),
        }
        for anc, acc in metrics["per_ancestry"].items():
            row[f"concordance_{anc}"] = round(acc, 4)

        rows.append(row)
        print(f"  Concordance {scenario}: {metrics['overall']:.2%}")

    return pd.DataFrame(rows)


def evaluate_global_ancestry(chrom: int) -> pd.DataFrame:
    """Compute mean global ancestry fractions for both tools on real data.

    Returns DataFrame with tool, scenario, and per-ancestry fractions.
    """
    rows = []

    for scenario in config.SCENARIOS:
        for tool_name in ["rfmix", "gnomix"]:
            if tool_name == "rfmix":
                msp_path = get_rfmix_msp_path(scenario, None, chrom)
            else:
                msp_path = get_gnomix_msp_path(scenario, None)

            if not msp_path.exists():
                continue

            result = parse_msp(msp_path, tool_name)
            intervals = msp_to_intervals(result)
            pred_df = intervals_to_df(intervals)

            fracs = mean_global_fractions(pred_df)

            row = {"tool": tool_name, "scenario": scenario}
            for anc, frac in fracs.items():
                row[f"frac_{anc}"] = round(frac, 4)

            rows.append(row)
            frac_str = ", ".join(f"{a}={f:.1%}" for a, f in fracs.items())
            print(f"  {tool_name} {scenario} global: {frac_str}")

    return pd.DataFrame(rows)


def run_evaluate(chrom: int) -> None:
    """Run all evaluation metrics and save to CSVs."""
    config.METRICS_DIR.mkdir(parents=True, exist_ok=True)

    print("Evaluating simulated data accuracy...")
    sim_df = evaluate_simulated(chrom)
    sim_path = config.METRICS_DIR / "simulated_accuracy.csv"
    sim_df.to_csv(sim_path, index=False)
    print(f"  Saved to {sim_path}")

    print("\nEvaluating real data concordance...")
    conc_df = evaluate_concordance(chrom)
    conc_path = config.METRICS_DIR / "concordance.csv"
    conc_df.to_csv(conc_path, index=False)
    print(f"  Saved to {conc_path}")

    print("\nComputing global ancestry fractions...")
    global_df = evaluate_global_ancestry(chrom)
    global_path = config.METRICS_DIR / "global_ancestry.csv"
    global_df.to_csv(global_path, index=False)
    print(f"  Saved to {global_path}")

    print("\nEvaluation complete.")
