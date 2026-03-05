"""Simulate admixed genomes with known ground-truth local ancestry.

Uses msprime for coalescent simulation with admixture demography,
and tspop to extract per-haplotype ancestry tracts as ground truth.
"""

import subprocess
import tempfile
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import tspop

import config


# Effective population sizes (standard values)
_NE = {"AFR": 14_474, "EUR": 10_000, "NAM": 10_000, "ADMIX": 30_000, "ANC": 14_474}

# Ancient divergence time (generations) for ancestral population
_T_DIVERGE = 2000

# Number of admixed diploid individuals to simulate
_N_ADMIXED = 50


def _build_demography(
    pop_labels: list[str],
    fractions: list[float],
    t_admix: int,
) -> tuple[msprime.Demography, float]:
    """Build an admixture demography for msprime.

    Creates source populations + an ADMIX population that forms at t_admix
    generations ago via mass migration from sources. A census event is placed
    at t_admix + 0.01 (required by tspop to identify ancestral lineages).
    An ancestral population (ANC) is added so all lineages can coalesce.

    Returns (demography, census_time).
    """
    demog = msprime.Demography()

    # Add admixed population (sampled at present)
    demog.add_population(name="ADMIX", initial_size=_NE["ADMIX"])

    # Add source populations
    for label in pop_labels:
        demog.add_population(name=label, initial_size=_NE[label])

    # Add ancestral population for coalescence
    demog.add_population(name="ANC", initial_size=_NE["ANC"])

    # Admixture events: use conditional proportions so all lineages leave ADMIX.
    # msprime applies mass migrations sequentially at same time:
    # after moving fraction[0], fraction[1] applies to the remainder, etc.
    # Last migration gets proportion=1.0 to drain all remaining lineages.
    remaining = 1.0
    for i, (label, frac) in enumerate(zip(pop_labels, fractions)):
        if i == len(pop_labels) - 1:
            prop = 1.0  # move all remaining
        else:
            prop = frac / remaining
            remaining -= frac
        demog.add_mass_migration(
            time=t_admix,
            source="ADMIX",
            dest=label,
            proportion=prop,
        )

    # Census event just above admixture time (needed by tspop)
    census_time = t_admix + 0.01
    demog.add_census(time=census_time)

    # Ancient divergence: all source populations merge into ANC
    for label in pop_labels:
        demog.add_mass_migration(
            time=_T_DIVERGE,
            source=label,
            dest="ANC",
            proportion=1.0,
        )

    demog.sort_events()
    return demog, census_time


def _load_recombination_map(chrom: int) -> msprime.RateMap:
    """Load an msprime RateMap from the HapMap-format genetic map."""
    return msprime.RateMap.read_hapmap(str(config.raw_genetic_map(chrom)))


def simulate_one(
    scenario: str,
    generation: int,
    chrom: int,
    seed: int | None = None,
) -> tuple[Path, Path]:
    """Simulate admixed individuals for one scenario x generation condition.

    Returns (vcf_path, ground_truth_path).
    """
    info = config.SCENARIOS[scenario]
    pop_labels = info["pops"]
    fractions = info["fractions"]

    demog, census_time = _build_demography(pop_labels, fractions, generation)
    rate_map = _load_recombination_map(chrom)

    rng = np.random.default_rng(seed)
    sim_seed = int(rng.integers(1, 2**31))
    mut_seed = int(rng.integers(1, 2**31))

    # Simulate ancestry
    ts = msprime.sim_ancestry(
        samples={"ADMIX": _N_ADMIXED},
        demography=demog,
        recombination_rate=rate_map,
        sequence_length=rate_map.sequence_length,
        random_seed=sim_seed,
    )

    # Add mutations
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=mut_seed)

    # Extract ground-truth ancestry tracts
    pop_anc = tspop.get_pop_ancestry(ts, census_time)
    gt_df = _extract_ground_truth(ts, pop_anc, pop_labels)

    # Write outputs
    out_dir = config.SIMULATED_DIR / scenario / f"gen_{generation}"
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = out_dir / "query.vcf.gz"
    gt_path = out_dir / "ground_truth.csv"

    _write_vcf(ts, chrom, vcf_path)
    gt_df.to_csv(gt_path, index=False)

    print(
        f"  Simulated {scenario}/gen_{generation}: {len(gt_df)} tracts, "
        f"{ts.num_mutations} variants"
    )

    return vcf_path, gt_path


def _extract_ground_truth(
    ts,
    pop_anc: tspop.PopAncestry,
    pop_labels: list[str],
) -> pd.DataFrame:
    """Convert tspop ancestry tracts to a ground-truth DataFrame.

    Maps population integer IDs back to string labels (AFR, EUR, NAM).
    Maps sample node IDs to individual/haplotype identifiers matching the VCF.

    Returns DataFrame with columns: [individual, haplotype, start, end, ancestry].
    """
    # Build pop ID -> label mapping from the tree sequence
    pop_id_to_label = {}
    for pop in ts.populations():
        name = pop.metadata.get("name", "") if isinstance(pop.metadata, dict) else ""
        if not name:
            name = str(pop.id)
        if name in pop_labels:
            pop_id_to_label[pop.id] = name

    squashed = pop_anc.squashed_table

    rows = []
    for _, tract in squashed.iterrows():
        node_id = int(tract["sample"])
        # Individual i has nodes 2*i (hap 0) and 2*i+1 (hap 1)
        indiv_idx = node_id // 2
        hap = node_id % 2
        indiv_name = f"tsk_{indiv_idx}"

        pop_id = int(tract["population"])
        label = pop_id_to_label.get(pop_id, f"UNK_{pop_id}")

        rows.append(
            {
                "individual": indiv_name,
                "haplotype": hap,
                "start": int(tract["left"]),
                "end": int(tract["right"]),
                "ancestry": label,
            }
        )

    return pd.DataFrame(rows)


def _write_vcf(ts, chrom: int, out_path: Path) -> None:
    """Write tree sequence to bgzipped, indexed VCF."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as tmp:
        ts.write_vcf(tmp, contig_id=str(chrom))
        tmp_path = tmp.name

    subprocess.run(["bgzip", "-f", tmp_path], check=True)
    bgz_path = tmp_path + ".gz"
    Path(bgz_path).rename(out_path)
    subprocess.run(["tabix", "-p", "vcf", str(out_path)], check=True)


def run_simulate(
    chrom: int,
    scenario: str | None = None,
    generations: list[int] | None = None,
    seed: int = 42,
) -> None:
    """Run simulation for specified or all scenario x generation conditions."""
    scenarios = [scenario] if scenario else list(config.SCENARIOS.keys())
    gens = generations if generations else config.ADMIXTURE_GENERATIONS

    rng = np.random.default_rng(seed)

    for sc in scenarios:
        for gen in gens:
            sim_seed = int(rng.integers(1, 2**31))
            print(f"Simulating {sc}/gen_{gen}...")
            simulate_one(sc, gen, chrom, seed=sim_seed)

    print("Simulation complete.")
