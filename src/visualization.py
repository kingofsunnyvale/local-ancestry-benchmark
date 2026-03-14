"""Visualization: all benchmark plots."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
from src.parsers import (
    get_rfmix_msp_path,
    get_gnomix_msp_path,
    intervals_to_df,
    msp_to_intervals,
    parse_msp,
)

# Consistent color palette for ancestries
ANC_COLORS = {"AFR": "#1b9e77", "EUR": "#d95f02", "NAM": "#7570b3"}
TOOL_COLORS = {"rfmix": "#e41a1c", "gnomix": "#377eb8"}
TOOL_LABELS = {"rfmix": "RFMix v2", "gnomix": "G-Nomix"}


def _savefig(fig, name: str) -> None:
    """Save figure to results/figures/ as both PDF and PNG."""
    config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "png"]:
        path = config.FIGURES_DIR / f"{name}.{ext}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {name}")


def plot_accuracy_vs_generations() -> None:
    """Line plot: overall accuracy vs admixture generations, one line per tool."""
    df = pd.read_csv(config.METRICS_DIR / "simulated_accuracy.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, scenario in zip(axes, ["2way", "3way"]):
        sub = df[df["scenario"] == scenario]
        for tool in ["rfmix", "gnomix"]:
            tsub = sub[sub["tool"] == tool].sort_values("generation")
            ax.plot(
                tsub["generation"],
                tsub["overall_accuracy"] * 100,
                marker="o",
                color=TOOL_COLORS[tool],
                label=TOOL_LABELS[tool],
                linewidth=2,
            )
        ax.set_xlabel("Admixture Generations Ago")
        ax.set_ylabel("Per-Base Accuracy (%)")
        ax.set_title(f"{scenario.upper()} Admixture")
        ax.legend()
        ax.set_xticks([5, 10, 20, 50])
        ax.grid(alpha=0.3)

    fig.suptitle("LAI Accuracy on Simulated Data", fontsize=14, y=1.02)
    _savefig(fig, "accuracy_vs_generations")


def plot_per_ancestry_accuracy() -> None:
    """Grouped bar chart: per-ancestry accuracy for each tool and scenario."""
    df = pd.read_csv(config.METRICS_DIR / "simulated_accuracy.csv")

    scenarios = ["2way", "3way"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, scenario in zip(axes, scenarios):
        sub = df[df["scenario"] == scenario]
        anc_cols = [c for c in sub.columns if c.startswith("accuracy_")]
        anc_names = [c.replace("accuracy_", "") for c in anc_cols]
        anc_names = [a for a in anc_names if sub[f"accuracy_{a}"].notna().any()]

        # Average across generations for each tool
        tools = ["rfmix", "gnomix"]
        x = np.arange(len(anc_names))
        width = 0.35

        for i, tool in enumerate(tools):
            tsub = sub[sub["tool"] == tool]
            means = [tsub[f"accuracy_{a}"].mean() * 100 for a in anc_names]
            ax.bar(
                x + i * width,
                means,
                width,
                label=TOOL_LABELS[tool],
                color=TOOL_COLORS[tool],
                alpha=0.85,
            )

        ax.set_xlabel("Ancestry")
        ax.set_ylabel("Mean Per-Base Accuracy (%)")
        ax.set_title(f"{scenario.upper()} — Per-Ancestry Accuracy")
        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(anc_names)
        ax.legend()
        ax.grid(alpha=0.3, axis="y")

    _savefig(fig, "per_ancestry_accuracy")


def plot_concordance() -> None:
    """Bar chart: concordance between tools on real data by scenario."""
    df = pd.read_csv(config.METRICS_DIR / "concordance.csv")

    fig, ax = plt.subplots(figsize=(6, 5))

    scenarios = df["scenario"].tolist()
    conc = df["concordance"].values * 100
    bars = ax.bar(scenarios, conc, color=["#66c2a5", "#fc8d62"], alpha=0.85, width=0.5)

    for bar, val in zip(bars, conc):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            f"{val:.1f}%",
            ha="center",
            va="bottom",
            fontsize=11,
        )

    ax.set_ylabel("Concordance (%)")
    ax.set_title("RFMix vs G-Nomix Concordance on Real Admixed Data")
    ax.set_ylim(0, 105)
    ax.grid(alpha=0.3, axis="y")
    _savefig(fig, "concordance")


def plot_global_ancestry() -> None:
    """Stacked bar chart: global ancestry fractions from both tools."""
    df = pd.read_csv(config.METRICS_DIR / "global_ancestry.csv")

    scenarios = df["scenario"].unique()
    fig, axes = plt.subplots(1, len(scenarios), figsize=(6 * len(scenarios), 5))
    if len(scenarios) == 1:
        axes = [axes]

    for ax, scenario in zip(axes, scenarios):
        sub = df[df["scenario"] == scenario]
        frac_cols = [c for c in sub.columns if c.startswith("frac_")]
        anc_names = [c.replace("frac_", "") for c in frac_cols]

        tools = sub["tool"].tolist()
        labels = [TOOL_LABELS.get(t, t) for t in tools]
        x = np.arange(len(labels))

        bottom = np.zeros(len(labels))
        for anc in anc_names:
            vals = sub[f"frac_{anc}"].fillna(0).values * 100
            ax.bar(
                x,
                vals,
                bottom=bottom,
                label=anc,
                color=ANC_COLORS.get(anc, "#999999"),
                width=0.5,
            )
            bottom += vals

        ax.set_ylabel("Global Ancestry (%)")
        ax.set_title(f"{scenario.upper()} — Global Ancestry Fractions")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(title="Ancestry")
        ax.set_ylim(0, 105)
        ax.grid(alpha=0.3, axis="y")

    _savefig(fig, "global_ancestry")


def plot_performance() -> None:
    """Grouped bar chart: wall-clock time per tool on simulated data."""
    df = pd.read_csv(config.METRICS_DIR / "performance.csv")
    sim = df[df["mode"] == "sim"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, scenario in zip(axes, ["2way", "3way"]):
        sub = sim[sim["scenario"] == scenario]
        tools = ["rfmix", "gnomix"]
        gens = sorted(sub["generation"].dropna().unique().astype(int))
        x = np.arange(len(gens))
        width = 0.35

        for i, tool in enumerate(tools):
            tsub = sub[sub["tool"] == tool].sort_values("generation")
            ax.bar(
                x + i * width,
                tsub["wall_seconds"].values,
                width,
                label=TOOL_LABELS[tool],
                color=TOOL_COLORS[tool],
                alpha=0.85,
            )

        ax.set_xlabel("Admixture Generations")
        ax.set_ylabel("Wall-Clock Time (s)")
        ax.set_title(f"{scenario.upper()} — Runtime")
        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(gens)
        ax.legend()
        ax.grid(alpha=0.3, axis="y")

    fig.suptitle("Computational Performance (Simulated Data)", fontsize=14, y=1.02)
    _savefig(fig, "performance_simulated")

    # Also plot real data times
    real = df[df["mode"] == "real"]
    if not real.empty:
        fig2, ax2 = plt.subplots(figsize=(6, 5))
        scenarios = real["scenario"].unique()
        x = np.arange(len(scenarios))
        width = 0.35

        for i, tool in enumerate(["rfmix", "gnomix"]):
            tsub = real[real["tool"] == tool]
            vals = [
                tsub[tsub["scenario"] == s]["wall_seconds"].values[0] for s in scenarios
            ]
            ax2.bar(
                x + i * width,
                vals,
                width,
                label=TOOL_LABELS[tool],
                color=TOOL_COLORS[tool],
                alpha=0.85,
            )

        ax2.set_xlabel("Scenario")
        ax2.set_ylabel("Wall-Clock Time (s)")
        ax2.set_title("Runtime on Real Admixed Data")
        ax2.set_xticks(x + width / 2)
        ax2.set_xticklabels([s.upper() for s in scenarios])
        ax2.legend()
        ax2.grid(alpha=0.3, axis="y")
        _savefig(fig2, "performance_real")


def plot_chromosome_painting(
    scenario: str = "2way",
    generation: int = 10,
    chrom: int = 22,
    n_individuals: int = 5,
) -> None:
    """Horizontal stacked bar showing ancestry along the chromosome for each haplotype."""
    fig, axes = plt.subplots(
        n_individuals * 2,
        2,
        figsize=(16, n_individuals * 2.5),
        sharex=True,
    )

    tools = ["rfmix", "gnomix"]

    for col, tool in enumerate(tools):
        if tool == "rfmix":
            msp_path = get_rfmix_msp_path(scenario, generation, chrom)
        else:
            msp_path = get_gnomix_msp_path(scenario, generation)

        if not msp_path.exists():
            continue

        result = parse_msp(msp_path, tool)
        intervals = msp_to_intervals(result)
        pred_df = intervals_to_df(intervals)

        individuals = sorted(pred_df["individual"].unique())[:n_individuals]

        row_idx = 0
        for indiv in individuals:
            for hap in [0, 1]:
                ax = axes[row_idx, col]
                hap_data = pred_df[
                    (pred_df["individual"] == indiv) & (pred_df["haplotype"] == hap)
                ]

                for _, iv in hap_data.iterrows():
                    color = ANC_COLORS.get(iv["ancestry"], "#999999")
                    ax.barh(
                        0,
                        iv["end"] - iv["start"],
                        left=iv["start"],
                        color=color,
                        height=0.8,
                        edgecolor="none",
                    )

                ax.set_yticks([])
                if col == 0:
                    ax.set_ylabel(
                        f"{indiv}.{hap}",
                        fontsize=8,
                        rotation=0,
                        ha="right",
                        va="center",
                    )
                ax.set_xlim(pred_df["start"].min(), pred_df["end"].max())

                if row_idx == 0:
                    ax.set_title(TOOL_LABELS[tool], fontsize=12)

                row_idx += 1

    # Add a legend
    from matplotlib.patches import Patch

    pops = sorted(ANC_COLORS.keys())
    legend_elements = [
        Patch(facecolor=ANC_COLORS[p], label=p)
        for p in pops
        if p in pred_df["ancestry"].values
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=len(pops),
        fontsize=10,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(
        f"Chromosome Painting — {scenario.upper()} gen={generation}",
        fontsize=14,
        y=1.01,
    )
    plt.tight_layout()
    _savefig(fig, f"chromosome_painting_{scenario}_gen{generation}")


def run_visualize(chrom: int = 22) -> None:
    """Generate all plots."""
    print("Generating visualizations...")

    plot_accuracy_vs_generations()
    plot_per_ancestry_accuracy()
    plot_concordance()
    plot_global_ancestry()
    plot_performance()

    # Chromosome paintings for a few conditions
    for scenario in ["2way", "3way"]:
        plot_chromosome_painting(scenario=scenario, generation=10, chrom=chrom)

    print(f"\nAll figures saved to {config.FIGURES_DIR}")
