"""Subprocess wrappers for RFMix and Gnomix with timing and memory tracking."""

import resource
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import config


@dataclass
class RunResult:
    """Stores tool execution metadata."""

    tool: str
    scenario: str
    generation: int | None  # None for real-data runs
    mode: str  # "sim" or "real"
    wall_seconds: float
    peak_memory_mb: float
    output_dir: Path


def _get_peak_memory_mb() -> float:
    """Get peak memory of child processes in MB (macOS returns bytes, Linux KB)."""
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    import sys

    if sys.platform == "darwin":
        return usage.ru_maxrss / (1024 * 1024)
    return usage.ru_maxrss / 1024


def run_rfmix(
    query_vcf: Path,
    reference_vcf: Path,
    sample_map: Path,
    genetic_map: Path,
    output_basename: Path,
    chrom: int,
) -> float:
    """Run RFMix and return wall-clock seconds."""
    output_basename.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(config.RFMIX_BIN),
        "-f",
        str(query_vcf),
        "-r",
        str(reference_vcf),
        "-m",
        str(sample_map),
        "-g",
        str(genetic_map),
        "-o",
        str(output_basename),
        f"--chromosome={chrom}",
    ]

    t0 = time.perf_counter()
    subprocess.run(cmd, check=True)
    return time.perf_counter() - t0


def run_gnomix(
    query_vcf: Path,
    reference_vcf: Path,
    sample_map: Path,
    genetic_map: Path,
    output_dir: Path,
    chrom: int,
) -> float:
    """Run Gnomix in training mode and return wall-clock seconds."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(config.GNOMIX_PYTHON),
        str(config.GNOMIX_SCRIPT),
        str(query_vcf),
        str(output_dir),
        str(chrom),
        "False",  # phase=False (data already phased)
        str(genetic_map),
        str(reference_vcf),
        str(sample_map),
    ]

    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, cwd=config.TOOLS_DIR / "gnomix")
    return time.perf_counter() - t0


def _run_tool_on_condition(
    tool: str,
    query_vcf: Path,
    scenario: str,
    generation: int | None,
    mode: str,
    chrom: int,
) -> RunResult:
    """Run a single tool on one query VCF with the matching reference data."""
    reference_vcf = config.PROCESSED_VCF_DIR / f"reference_{scenario}_chr{chrom}.vcf.gz"
    sample_map = config.PROCESSED_SMAPS_DIR / f"reference_{scenario}.txt"
    genetic_map = config.PROCESSED_GMAPS_DIR / f"chr{chrom}.txt"

    if mode == "sim":
        tag = f"{scenario}/gen_{generation}"
    else:
        tag = f"{scenario}/real"

    mem_before = _get_peak_memory_mb()

    if tool == "rfmix":
        out_dir = config.RFMIX_DIR / tag
        out_basename = out_dir / f"chr{chrom}"
        elapsed = run_rfmix(
            query_vcf,
            reference_vcf,
            sample_map,
            genetic_map,
            out_basename,
            chrom,
        )
    elif tool == "gnomix":
        out_dir = config.GNOMIX_DIR / tag
        elapsed = run_gnomix(
            query_vcf,
            reference_vcf,
            sample_map,
            genetic_map,
            out_dir,
            chrom,
        )
    else:
        raise ValueError(f"Unknown tool: {tool}")

    mem_after = _get_peak_memory_mb()
    peak_mb = max(0, mem_after - mem_before)

    result = RunResult(
        tool=tool,
        scenario=scenario,
        generation=generation,
        mode=mode,
        wall_seconds=elapsed,
        peak_memory_mb=peak_mb,
        output_dir=out_dir,
    )

    print(f"  {tool} on {tag}: {elapsed:.1f}s, ~{peak_mb:.0f} MB peak memory")
    return result


def run_tools(
    chrom: int,
    tool: str | None = None,
    mode: str | None = None,
    scenario: str | None = None,
    generations: list[int] | None = None,
) -> list[RunResult]:
    """Run LAI tools on simulated and/or real query data.

    Parameters
    ----------
    chrom : Chromosome number.
    tool : "rfmix", "gnomix", or None (both).
    mode : "sim", "real", or None (both).
    scenario : "2way", "3way", or None (all).
    generations : List of generation values, or None (all from config).
    """
    tools = [tool] if tool else ["rfmix", "gnomix"]
    scenarios = [scenario] if scenario else list(config.SCENARIOS.keys())
    gens = generations if generations else config.ADMIXTURE_GENERATIONS
    modes = [mode] if mode else ["sim", "real"]

    results = []

    for t in tools:
        for sc in scenarios:
            if "sim" in modes:
                for gen in gens:
                    query_vcf = (
                        config.SIMULATED_DIR / sc / f"gen_{gen}" / "query.vcf.gz"
                    )
                    if not query_vcf.exists():
                        print(f"  Skipping {t} on {sc}/gen_{gen}: query VCF not found")
                        continue
                    print(f"Running {t} on {sc}/gen_{gen}...")
                    r = _run_tool_on_condition(t, query_vcf, sc, gen, "sim", chrom)
                    results.append(r)

            if "real" in modes:
                query_vcf = config.PROCESSED_VCF_DIR / f"query_real_chr{chrom}.vcf.gz"
                if not query_vcf.exists():
                    print(f"  Skipping {t} on {sc}/real: query VCF not found")
                    continue
                print(f"Running {t} on {sc}/real...")
                r = _run_tool_on_condition(t, query_vcf, sc, None, "real", chrom)
                results.append(r)

    # Save performance summary
    _save_performance(results)
    print(f"Tool execution complete. {len(results)} runs finished.")
    return results


def _save_performance(results: list[RunResult]) -> None:
    """Save timing and memory results to CSV."""
    if not results:
        return
    config.METRICS_DIR.mkdir(parents=True, exist_ok=True)
    import csv

    out = config.METRICS_DIR / "performance.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "tool",
                "scenario",
                "generation",
                "mode",
                "wall_seconds",
                "peak_memory_mb",
            ],
        )
        writer.writeheader()
        for r in results:
            writer.writerow(
                {
                    "tool": r.tool,
                    "scenario": r.scenario,
                    "generation": r.generation if r.generation else "",
                    "mode": r.mode,
                    "wall_seconds": round(r.wall_seconds, 2),
                    "peak_memory_mb": round(r.peak_memory_mb, 1),
                }
            )
    print(f"  Performance metrics saved to {out}")
