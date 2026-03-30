# ══════════════════════════════════════════════════════════════
# Snakefile — AD RNA-seq Analysis Pipeline
# Dataset: GSE125583 — Alzheimer's Disease Fusiform Gyrus
# ══════════════════════════════════════════════════════════════
#
# Usage:
#   snakemake -n                          # dry run
#   snakemake --cores 4                   # run locally
#   snakemake --cores 4 --forceall        # force re-run everything
#   snakemake --dag | dot -Tpng > dag.png # visualize pipeline
#
# SLURM usage:
#   sbatch --wrap="module load python/3.13.5 && \
#     source ~/rnaseq_venv311/bin/activate && \
#     cd /**project path** && \
#     snakemake --cores 4" \
#     --output=logs/snakemake_%j.out --error=logs/snakemake_%j.err \
#     --time=04:00:00 --mem=32G --cpus-per-task=4
# ══════════════════════════════════════════════════════════════

import glob

# ── Configuration ─────────────────────────────────────────────
RAW_DIR      = "data/raw"
META_DIR     = "data/metadata"
PROC_DIR     = "data/processed"
ANNOT_DIR    = "data/annotations"
RESULTS_DIR  = "analysis_v2"
QC_DIR       = "results/qc"

SOFT_FILE    = f"{RAW_DIR}/GSE125583_family.soft.gz"
GENE_INFO    = f"{ANNOT_DIR}/Homo_sapiens.gene_info.gz"

# ── Target rule ───────────────────────────────────────────────
rule all:
    input:
        # Script 1 outputs
        f"{RAW_DIR}/merged_counts.tsv",
        f"{META_DIR}/sample_metadata_full.csv",
        # Script 2 outputs
        f"{PROC_DIR}/02_filtered_counts.csv",
        f"{QC_DIR}/01_library_sizes.pdf",
        # Script 3 outputs
        f"{ANNOT_DIR}/entrez_to_symbol.csv",
        f"{PROC_DIR}/04_deseq2_normalized_counts_symbols.csv",
        f"{META_DIR}/sample_metadata_covariates.csv",
        # Script 4 outputs
        f"{RESULTS_DIR}/results_model1/DE_results.csv",
        # Script 5 outputs
        f"{RESULTS_DIR}/results_model2/DE_results.csv",
        f"{RESULTS_DIR}/results_model3/DE_results.csv",
        f"{RESULTS_DIR}/results_model4/DE_results.csv",
        # Script 6 outputs
        f"{RESULTS_DIR}/verification_report.txt",
        f"{RESULTS_DIR}/enrichment/BIOLOGICAL_INTERPRETATION.txt",

# ── Rule 1: Merge counts + build metadata ────────────────────
rule merge_and_build_metadata:
    input:
        soft = SOFT_FILE,
        samples = glob.glob(f"{RAW_DIR}/GSM*.tsv.gz"),
    output:
        counts = f"{RAW_DIR}/merged_counts.tsv",
        meta   = f"{META_DIR}/sample_metadata_full.csv",
    log:
        "logs/01_merge.log"
    shell:
        "python scripts/01_merge_and_build_metadata.py > {log} 2>&1"

# ── Rule 2: QC and filter ────────────────────────────────────
rule qc_and_filter:
    input:
        counts = f"{RAW_DIR}/merged_counts.tsv",
        meta   = f"{META_DIR}/sample_metadata_full.csv",
    output:
        filtered = f"{PROC_DIR}/02_filtered_counts.csv",
        filt_meta = f"{PROC_DIR}/02_filtered_metadata.csv",
        qc_plot  = f"{QC_DIR}/01_library_sizes.pdf",
    log:
        "logs/02_qc_filter.log"
    shell:
        "python scripts/02_qc_and_filter.py > {log} 2>&1"

# ── Rule 3: Build covariates ─────────────────────────────────
rule build_covariates:
    input:
        filtered  = f"{PROC_DIR}/02_filtered_counts.csv",
        meta      = f"{META_DIR}/sample_metadata_full.csv",
        gene_info = GENE_INFO,
    output:
        mapping    = f"{ANNOT_DIR}/entrez_to_symbol.csv",
        norm_sym   = f"{PROC_DIR}/04_deseq2_normalized_counts_symbols.csv",
        norm_ent   = f"{PROC_DIR}/04_deseq2_normalized_counts.csv",
        covariates = f"{META_DIR}/sample_metadata_covariates.csv",
    log:
        "logs/03_covariates.log"
    shell:
        "python scripts/03_build_covariates.py > {log} 2>&1"

# ── Rule 4: Baseline DESeq2 (Model 1) ────────────────────────
rule baseline_deseq2:
    input:
        counts = f"{PROC_DIR}/02_filtered_counts.csv",
        meta   = f"{META_DIR}/sample_metadata_covariates.csv",
    output:
        results = f"{RESULTS_DIR}/results_model1/DE_results.csv",
    log:
        "logs/04_baseline_deseq2.log"
    threads: 4
    shell:
        "python scripts/04_baseline_deseq2.py > {log} 2>&1"

# ── Rule 5: All models DESeq2 (Models 2-4) + symbol fix ──────
rule all_models_deseq2:
    input:
        counts  = f"{PROC_DIR}/02_filtered_counts.csv",
        meta    = f"{META_DIR}/sample_metadata_covariates.csv",
        mapping = f"{ANNOT_DIR}/entrez_to_symbol.csv",
        model1  = f"{RESULTS_DIR}/results_model1/DE_results.csv",
    output:
        model2 = f"{RESULTS_DIR}/results_model2/DE_results.csv",
        model3 = f"{RESULTS_DIR}/results_model3/DE_results.csv",
        model4 = f"{RESULTS_DIR}/results_model4/DE_results.csv",
    log:
        "logs/05_all_models.log"
    threads: 4
    shell:
        "python scripts/05_all_models_deseq2.py > {log} 2>&1"

# ── Rule 6: Enrichment + verification ────────────────────────
rule enrichment_and_verification:
    input:
        model4 = f"{RESULTS_DIR}/results_model4/DE_results.csv",
        model1 = f"{RESULTS_DIR}/results_model1/DE_results.csv",
        model2 = f"{RESULTS_DIR}/results_model2/DE_results.csv",
        model3 = f"{RESULTS_DIR}/results_model3/DE_results.csv",
        meta   = f"{META_DIR}/sample_metadata_covariates.csv",
    output:
        report  = f"{RESULTS_DIR}/verification_report.txt",
        interp  = f"{RESULTS_DIR}/enrichment/BIOLOGICAL_INTERPRETATION.txt",
    log:
        "logs/06_enrichment.log"
    threads: 4
    shell:
        "python scripts/06_enrichment_and_verification.py > {log} 2>&1"

# ── Utility: clean all generated files ────────────────────────
rule clean:
    shell:
        """
        rm -rf {PROC_DIR}/*
        rm -rf {META_DIR}/sample_metadata_full.csv
        rm -rf {META_DIR}/sample_metadata_covariates.csv
        rm -rf {ANNOT_DIR}/entrez_to_symbol.csv
        rm -rf {RAW_DIR}/merged_counts.tsv
        rm -rf {RESULTS_DIR}/results_model*/
        rm -rf {RESULTS_DIR}/enrichment/
        rm -rf {RESULTS_DIR}/verification_report.txt
        rm -rf {QC_DIR}/*
        rm -rf logs/*.log
        echo "Cleaned all generated files."
        """
