# Local Ancestry Inference Benchmark

Benchmarking pipeline comparing two local ancestry inference methods, RFMIX and GNOMIX, on admixed human genomes using chromosome 22 (for now)from the 1000 Genomes Project.

## Background

Local ancestry inference assigns ancestry labels to segments of an individual's genome, identifying which ancestral population each chromosomal region was inherited from. This is important for admixture mapping, correcting for population stratification in GWAS, and understanding demographic history.

We look at two admixture scenarios:

- **2-way admixture** (AFR + EUR) — modeling African American ancestry
- **3-way admixture** (AFR + EUR + NAM) — modeling Latino ancestry

## Reference and Query Populations

| Role | Superpopulation | 1KG Subpopulation |
|------|----------------|-------------------|
| Reference | AFR | YRI (Yoruba) |
| Reference | EUR | CEU (Utah European) |
| Reference | NAM | PEL (Peruvian) |
| Query (real) | Admixed | ASW, ACB, MXL, PUR |

## Steps

1. **`prep`** — Reformat genetic maps, create sample maps, split VCFs into reference and query panels
2. **`simulate`** — Generate admixed genomes with known ground-truth ancestry tracts using msprime + tspop
3. **`run`** — Execute RFMIX and GNOMIX on simulated and real query data, recording runtime and memory
4. **`evaluate`** — Parse tool outputs, compute per-base accuracy (simulated), inter-tool concordance (real), and global ancestry fractions
5. **`visualize`** — Generate figures: accuracy curves, confusion matrices, concordance plots, chromosome paintings, performance comparisons
6. **`robustness`** — Test sensitivity to reference panel size by downsampling (25%, 50%, 75%, 100%)

## Usage

All commands are run via `uv run python main.py <command>`. Every subcommand accepts `--chrom` (default: 22).

### Run the full pipeline

```bash
uv run python main.py all --chrom 22
```

### Run individual steps

```bash
# 1. Prepare data (reformat genetic maps, create sample maps, split VCFs)
uv run python main.py prep --chrom 22

# 2. Simulate admixed genomes with ground-truth ancestry
uv run python main.py simulate --chrom 22
uv run python main.py simulate --chrom 22 --scenario 2way --generations 5,10

# 3. Run LAI tools on query data
uv run python main.py run --chrom 22
uv run python main.py run --chrom 22 --tool rfmix --mode sim
uv run python main.py run --chrom 22 --tool gnomix --mode real --scenario 3way

# 4. Evaluate results (accuracy, concordance, global ancestry fractions)
uv run python main.py evaluate --chrom 22

# 5. Generate all plots
uv run python main.py visualize --chrom 22
```

### Key options

| Option | Subcommands | Description |
|--------|------------|-------------|
| `--chrom` | all | Chromosome number (default: 22) |
| `--scenario` | simulate, run, all | `2way` or `3way` (default: both) |
| `--generations` | simulate, run, all | Comma-separated list, e.g. `5,10,20,50` (default: all) |
| `--tool` | run, all | `rfmix` or `gnomix` (default: both) |
| `--mode` | run, all | `sim` or `real` (default: both) |
| `--seed` | simulate, all | Random seed for simulation (default: 42) |