# Alzheimer's Disease Bulk RNA-seq Analysis

**Author:** Durga Gomathi Arumuganainar  
**Date:** March 2026  
**Dataset:** [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583) — Fusiform Gyrus, 289 samples (219 AD, 70 Control)  
**Reference:** Srinivasan et al. 2020, *Cell Reports* 31(13):107843

---

## Key Finding

Naive differential expression identifies 11,878 genes (59%). After correcting for age, sex, Braak stage, and cell-type composition, only 2,289 genes (11.4%) remain — demonstrating that the majority of bulk RNA-seq signal in AD brain reflects cellular composition changes, not transcriptional regulation.

| Model | Covariates | Significant Genes | % |
|-------|-----------|------------------:|-----:|
| 1 | Diagnosis only | 11,878 | 59.3 |
| 2 | + Age, Sex | 11,459 | 57.2 |
| 3 | + Braak stage | 7,202 | 36.0 |
| **4** | **+ Cell composition** | **2,289** | **11.4** |

---

## Pipeline

```
01_merge_and_build_metadata.py   → Merge 289 count files + parse GEO metadata
02_qc_and_filter.py              → Filter low-count genes (30,727 → 20,024)
03_build_covariates.py           → Entrez→Symbol mapping + cell-type scores
04_baseline_deseq2.py            → Model 1: ~ Diagnosis
05_all_models_deseq2.py          → Models 2–4 + gene symbol fix
06_enrichment_and_verification.py → GO/GSEA enrichment + biological validation
Snakefile                        → Reproducible pipeline wrapper
```

---

## Quick Start

```bash
# Environment
module load python/3.13.5
source ~/rnaseq_venv311/bin/activate

# Dry run
snakemake -n

# Full run (submit via SLURM)
sbatch --wrap="snakemake --cores 4" \
       --time=04:00:00 --mem=32G --cpus-per-task=4
```

---

## Project Structure

```
AD_RNAseq_clean/
├── scripts/                  # 6 analysis scripts
├── data/
│   ├── raw/                  # 289 GSM files + SOFT + merged counts
│   ├── metadata/             # Clinical metadata + covariates
│   ├── processed/            # Filtered + normalized counts
│   └── annotations/          # Entrez→Symbol mapping
├── analysis_v2/
│   ├── results_model1–4/     # DESeq2 results per model
│   ├── enrichment/           # GO, GSEA, KEGG results
│   └── verification_report.txt
├── results/qc/               # QC plots
├── Snakefile
└── README.md
```

---

## Dependencies

`pandas` · `numpy` · `matplotlib` · `pydeseq2` · `gseapy` · `scipy` · `snakemake`
