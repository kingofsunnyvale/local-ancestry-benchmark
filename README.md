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