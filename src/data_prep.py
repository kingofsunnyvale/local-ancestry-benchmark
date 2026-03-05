"""Data preparation: genetic map reformatting, sample map creation, VCF splitting."""

import csv
import subprocess
from pathlib import Path

import config


def reformat_genetic_map(chrom: int) -> Path:
    """Convert 4-column genetic map to 3-column (Chromosome, Position(bp), Map(cM)).

    Strips 'chr' prefix from chromosome column. Both RFMix and Gnomix need
    this 3-column format; without it RFMix silently misreads columns and
    Gnomix crashes.
    """
    config.PROCESSED_GMAPS_DIR.mkdir(parents=True, exist_ok=True)
    src = config.raw_genetic_map(chrom)
    dst = config.PROCESSED_GMAPS_DIR / f"chr{chrom}.txt"

    with open(src) as fin, open(dst, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        next(reader)  # skip header
        fout.write("chr\tpos\tMap(cM)\n")
        for row in reader:
            chromosome = row[0].replace("chr", "")
            position = row[1]
            map_cm = row[3]
            fout.write(f"{chromosome}\t{position}\t{map_cm}\n")

    print(f"  Genetic map: {dst}")
    return dst


def _load_panel() -> list[dict]:
    """Read the 1KG panel file into a list of dicts."""
    rows = []
    with open(config.PANEL_FILE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def create_sample_maps() -> dict[str, Path]:
    """Create sample map files for 2-way and 3-way reference panels.

    Each sample map is a 2-column TSV: sample_id <tab> superpop_label.
    Labels are the superpop names (AFR, EUR, NAM) so both tools produce
    identical integer mappings (alphabetical order).

    Returns dict mapping scenario name to sample map path.
    """
    config.PROCESSED_SMAPS_DIR.mkdir(parents=True, exist_ok=True)
    panel = _load_panel()

    # Build lookup: subpop -> superpop label
    subpop_to_label = {}
    for label, subpops in config.REFERENCE_POPS.items():
        for sp in subpops:
            subpop_to_label[sp] = label

    results = {}
    for scenario, info in config.SCENARIOS.items():
        wanted_labels = set(info["pops"])
        dst = config.PROCESSED_SMAPS_DIR / f"reference_{scenario}.txt"

        with open(dst, "w") as f:
            for row in panel:
                pop = row["pop"]
                if pop in subpop_to_label:
                    label = subpop_to_label[pop]
                    if label in wanted_labels:
                        f.write(f"{row['sample']}\t{label}\n")

        results[scenario] = dst
        print(f"  Sample map ({scenario}): {dst}")

    return results


def _get_samples_for_pops(panel: list[dict], pops: list[str]) -> list[str]:
    """Get sample IDs belonging to the given subpopulation codes."""
    return [row["sample"] for row in panel if row["pop"] in pops]


def split_vcfs(chrom: int) -> dict[str, Path]:
    """Split the 1KG VCF into reference and query VCFs.

    Creates:
    - reference_2way_chr{chrom}.vcf.gz (YRI + CEU)
    - reference_3way_chr{chrom}.vcf.gz (YRI + CEU + PEL)
    - query_real_chr{chrom}.vcf.gz (ASW + ACB + MXL + PUR)

    Returns dict mapping name to path.
    """
    config.PROCESSED_VCF_DIR.mkdir(parents=True, exist_ok=True)
    panel = _load_panel()
    src_vcf = config.raw_vcf(chrom)

    outputs = {}

    # Reference VCFs per scenario
    for scenario, info in config.SCENARIOS.items():
        ref_subpops = []
        for label in info["pops"]:
            ref_subpops.extend(config.REFERENCE_POPS[label])
        samples = _get_samples_for_pops(panel, ref_subpops)
        name = f"reference_{scenario}_chr{chrom}"
        dst = config.PROCESSED_VCF_DIR / f"{name}.vcf.gz"
        _subset_vcf(src_vcf, samples, dst)
        outputs[f"reference_{scenario}"] = dst
        print(f"  VCF ({name}): {dst}")

    # Query VCF (real admixed individuals)
    query_samples = _get_samples_for_pops(panel, config.ADMIXED_POPS)
    name = f"query_real_chr{chrom}"
    dst = config.PROCESSED_VCF_DIR / f"{name}.vcf.gz"
    _subset_vcf(src_vcf, query_samples, dst)
    outputs["query_real"] = dst
    print(f"  VCF ({name}): {dst}")

    return outputs


def _subset_vcf(src: Path, samples: list[str], dst: Path) -> None:
    """Subset a VCF to the given samples using bcftools, then index."""
    sample_str = ",".join(samples)
    cmd_view = [
        "bcftools",
        "view",
        "-S",
        f"<(echo '{sample_str}' | tr ',' '\\n')",
        "--force-samples",
        str(src),
        "-Oz",
        "-o",
        str(dst),
    ]
    # bcftools -S expects a file; use process substitution via bash
    subprocess.run(
        ["bash", "-c", " ".join(cmd_view)],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(dst)], check=True)


def run_prep(chrom: int) -> None:
    """Run all data preparation steps."""
    print(f"Preparing data for chr{chrom}...")

    print("Reformatting genetic map...")
    reformat_genetic_map(chrom)

    print("Creating sample maps...")
    create_sample_maps()

    print("Splitting VCFs...")
    split_vcfs(chrom)

    print("Data preparation complete.")
