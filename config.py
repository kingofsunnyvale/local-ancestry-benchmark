from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

# Raw data
DATA_DIR = PROJECT_ROOT / "data"
GENETIC_MAPS_DIR = DATA_DIR / "genetic_maps"
PANELS_DIR = DATA_DIR / "panels"
VCF_DIR = DATA_DIR / "vcf"

# Processed data
PROCESSED_DIR = DATA_DIR / "processed"
PROCESSED_GMAPS_DIR = PROCESSED_DIR / "genetic_maps"
PROCESSED_SMAPS_DIR = PROCESSED_DIR / "sample_maps"
PROCESSED_VCF_DIR = PROCESSED_DIR / "vcf"

# Results
RESULTS_DIR = PROJECT_ROOT / "results"
SIMULATED_DIR = RESULTS_DIR / "simulated"
RFMIX_DIR = RESULTS_DIR / "rfmix"
GNOMIX_DIR = RESULTS_DIR / "gnomix"
METRICS_DIR = RESULTS_DIR / "metrics"
FIGURES_DIR = RESULTS_DIR / "figures"

# Tools
TOOLS_DIR = PROJECT_ROOT / "tools"
RFMIX_BIN = TOOLS_DIR / "rfmix" / "rfmix"
GNOMIX_SCRIPT = TOOLS_DIR / "gnomix" / "gnomix.py"
GNOMIX_PYTHON = TOOLS_DIR / "gnomix" / ".venv" / "bin" / "python"

# Panel file
PANEL_FILE = PANELS_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

# Reference populations: superpop label -> list of subpops to include
REFERENCE_POPS = {
    "AFR": ["YRI"],
    "EUR": ["CEU"],
    "NAM": ["PEL"],
}

# Admixed populations to use as real-data queries
ADMIXED_POPS = ["ASW", "ACB", "MXL", "PUR"]

# Simulation parameters
ADMIXTURE_GENERATIONS = [5, 10, 20, 50]

# Robustness: reference panel downsample fractions
DOWNSAMPLE_FRACTIONS = [0.25, 0.5, 0.75, 1.0]

# Scenarios
SCENARIOS = {
    "2way": {"pops": ["AFR", "EUR"], "fractions": [0.8, 0.2]},
    "3way": {"pops": ["AFR", "EUR", "NAM"], "fractions": [0.5, 0.3, 0.2]},
}


def raw_vcf(chrom: int) -> Path:
    return (
        VCF_DIR
        / f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    )


def raw_genetic_map(chrom: int) -> Path:
    return GENETIC_MAPS_DIR / f"genetic_map_GRCh37_chr{chrom}.txt"
