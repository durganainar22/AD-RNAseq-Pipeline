#!/usr/bin/env python
"""
04_baseline_deseq2.py
----------------------
Run baseline (naive) DESeq2: ~ Diagnosis only.
This is Model 1 — no covariates, just AD vs Control.
Expected result: ~11,000-12,000 significant genes (heavily confounded).

Inputs:
    data/processed/02_filtered_counts.csv           (20,024 genes x 289 samples)
    data/metadata/sample_metadata_covariates.csv    (289 samples with covariates)

Outputs:
    analysis_v2/results_model1/DE_results.csv       (all genes with stats)
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path('/scratch/arumuganainar.d/projects/AD_RNAseq_Project/AD_RNAseq_clean')
os.chdir(BASE)

print("=" * 65)
print("SCRIPT 4: Baseline DESeq2 (Model 1)")
print("  Design: ~ Diagnosis")
print("  No covariates — naive AD vs Control comparison")
print("=" * 65)

# ═════════════════════════════════════════════════════════════
# Load data
# ═════════════════════════════════════════════════════════════
print("\n--- Loading data ---\n")

counts = pd.read_csv(BASE / 'data' / 'processed' / '02_filtered_counts.csv', index_col=0)
meta = pd.read_csv(BASE / 'data' / 'metadata' / 'sample_metadata_covariates.csv')

# PyDESeq2 needs: samples as rows, genes as columns, integer counts
counts_int = counts.astype(int).T  # transpose: samples x genes

# Align metadata to counts
meta = meta.set_index('SampleID')
meta = meta.loc[counts_int.index].copy()

# Diagnosis must be categorical
meta['Diagnosis'] = meta['Diagnosis'].astype('category')

print(f"Count matrix: {counts_int.shape[0]} samples x {counts_int.shape[1]:,} genes")
print(f"\nDiagnosis distribution:")
print(meta['Diagnosis'].value_counts().to_string())

# ═════════════════════════════════════════════════════════════
# Run DESeq2 — Model 1: ~ Diagnosis
# ═════════════════════════════════════════════════════════════
print("\n--- Running DESeq2 (Model 1: ~ Diagnosis) ---")
print("This will take ~10-15 minutes...\n")

dds = DeseqDataSet(
    counts=counts_int,
    metadata=meta,
    design_factors=['Diagnosis'],
    refit_cooks=True,
)
dds.deseq2()

print("DESeq2 fitting complete. Running Wald test...\n")

stat_res = DeseqStats(dds, contrast=['Diagnosis', 'AD', 'Control'])
stat_res.summary()

# ═════════════════════════════════════════════════════════════
# Process and save results
# ═════════════════════════════════════════════════════════════
print("\n--- Processing results ---\n")

res = stat_res.results_df.copy()

# Add significance and direction columns
res['Significant'] = res['padj'] < 0.05
res['Direction'] = 'NS'
res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] > 0), 'Direction'] = 'Up'
res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] < 0), 'Direction'] = 'Down'

# Summary
n_total = len(res)
n_sig = res['Significant'].sum()
n_up = (res['Direction'] == 'Up').sum()
n_down = (res['Direction'] == 'Down').sum()
pct_sig = n_sig / n_total * 100

print(f"Total genes tested: {n_total:,}")
print(f"Significant (padj < 0.05): {n_sig:,} ({pct_sig:.1f}%)")
print(f"  Upregulated in AD:   {n_up:,}")
print(f"  Downregulated in AD: {n_down:,}")

# Top genes
print(f"\nTop 10 most significant genes:")
top = res.dropna(subset=['padj']).sort_values('padj').head(10)
for gene, row in top.iterrows():
    print(f"  {str(gene):<15} LFC={row['log2FoldChange']:7.3f}  "
          f"padj={row['padj']:.2e}  {row['Direction']}")

# Save
outdir = BASE / 'analysis_v2' / 'results_model1'
outdir.mkdir(parents=True, exist_ok=True)

out_file = outdir / 'DE_results.csv'
res.to_csv(out_file)
print(f"\nSaved: {out_file}")

print("\n" + "=" * 65)
print("DONE — Model 1 (Baseline)")
print(f"  {n_sig:,} significant genes ({pct_sig:.1f}%) — EXPECTED to be inflated")
print(f"  This reflects confounding from cell composition + disease severity")
print("  Next: python scripts/05_all_models_deseq2.py")
print("=" * 65)
