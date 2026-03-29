#!/usr/bin/env python
"""
05_all_models_deseq2.py
------------------------
Run Models 2-4 with progressive covariate adjustment, then fix gene
symbols on ALL model results (1-4).

Model 2: ~ Age_c + Sex_binary + Diagnosis
Model 3: ~ Age_c + Sex_binary + Braak_c + Diagnosis
Model 4: ~ Age_c + Sex_binary + Braak_c + MicrogliaScore + NeuronScore
           + AstrocyteScore + OligoScore + Diagnosis

Inputs:
    data/processed/02_filtered_counts.csv
    data/metadata/sample_metadata_covariates.csv
    data/annotations/entrez_to_symbol.csv
    analysis_v2/results_model1/DE_results.csv      (from Script 4)

Outputs:
    analysis_v2/results_model2/DE_results.csv
    analysis_v2/results_model3/DE_results.csv
    analysis_v2/results_model4/DE_results.csv
    (also re-saves model1 with gene symbols)
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
print("SCRIPT 5: All DESeq2 Models (2-4) + Gene Symbol Fix")
print("=" * 65)

# ═════════════════════════════════════════════════════════════
# Load data
# ═════════════════════════════════════════════════════════════
print("\n--- Loading data ---\n")

counts = pd.read_csv(BASE / 'data' / 'processed' / '02_filtered_counts.csv', index_col=0)
meta = pd.read_csv(BASE / 'data' / 'metadata' / 'sample_metadata_covariates.csv')
mapping = pd.read_csv(BASE / 'data' / 'annotations' / 'entrez_to_symbol.csv')

# PyDESeq2 needs: samples as rows, genes as columns, integer counts
counts_int = counts.astype(int).T

# Align metadata
meta = meta.set_index('SampleID')
meta = meta.loc[counts_int.index].copy()

# Ensure types
meta['Diagnosis'] = meta['Diagnosis'].astype('category')
for col in ['Age_c', 'Sex_binary', 'Braak_c',
            'MicrogliaScore', 'NeuronScore', 'AstrocyteScore', 'OligoScore']:
    if col in meta.columns:
        meta[col] = pd.to_numeric(meta[col], errors='coerce')

print(f"Counts: {counts_int.shape[0]} samples x {counts_int.shape[1]:,} genes")
print(f"Diagnosis: {dict(meta['Diagnosis'].value_counts())}")

# ═════════════════════════════════════════════════════════════
# Helper function to run a model
# ═════════════════════════════════════════════════════════════
def run_model(name, design_factors, outdir):
    """Run DESeq2 with given design, dropping NAs in required covariates."""
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'─'*65}")
    print(f"  {name}")
    print(f"  Design: ~ {' + '.join(design_factors)}")
    print(f"{'─'*65}")

    # Drop samples with NA in required covariates
    needed = [x for x in design_factors if x != 'Diagnosis']
    meta_sub = meta.dropna(subset=needed) if needed else meta.copy()
    counts_sub = counts_int.loc[meta_sub.index]

    print(f"  Samples: {len(meta_sub)} (dropped {len(meta) - len(meta_sub)} with missing covariates)")
    print(f"  AD: {(meta_sub['Diagnosis'] == 'AD').sum()}  Control: {(meta_sub['Diagnosis'] == 'Control').sum()}")

    # Fit
    dds = DeseqDataSet(
        counts=counts_sub,
        metadata=meta_sub,
        design_factors=design_factors,
        refit_cooks=True,
    )
    dds.deseq2()

    # Wald test
    stat_res = DeseqStats(dds, contrast=['Diagnosis', 'AD', 'Control'])
    stat_res.summary()

    # Process results
    res = stat_res.results_df.copy()
    res['Significant'] = res['padj'] < 0.05
    res['Direction'] = 'NS'
    res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] > 0), 'Direction'] = 'Up'
    res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] < 0), 'Direction'] = 'Down'

    # Save (Entrez IDs as index for now — symbols fixed later)
    res.to_csv(outdir / 'DE_results.csv')

    n_sig = res['Significant'].sum()
    n_total = len(res)
    n_up = (res['Direction'] == 'Up').sum()
    n_down = (res['Direction'] == 'Down').sum()

    print(f"\n  Results: {n_sig:,} significant ({n_sig/n_total*100:.1f}%)")
    print(f"  Up: {n_up:,}  Down: {n_down:,}")

    return n_sig, n_total

# ═════════════════════════════════════════════════════════════
# Run Models 2-4
# ═════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("RUNNING MODELS (each takes ~10-15 min)")
print("=" * 65)

m2_sig, total = run_model(
    "Model 2 — Age + Sex",
    ['Age_c', 'Sex_binary', 'Diagnosis'],
    BASE / 'analysis_v2' / 'results_model2',
)

m3_sig, _ = run_model(
    "Model 3 — Age + Sex + Braak",
    ['Age_c', 'Sex_binary', 'Braak_c', 'Diagnosis'],
    BASE / 'analysis_v2' / 'results_model3',
)

m4_sig, _ = run_model(
    "Model 4 — Age + Sex + Braak + Cell Composition (FINAL)",
    ['Age_c', 'Sex_binary', 'Braak_c',
     'MicrogliaScore', 'NeuronScore', 'AstrocyteScore', 'OligoScore',
     'Diagnosis'],
    BASE / 'analysis_v2' / 'results_model4',
)

# ═════════════════════════════════════════════════════════════
# Fix gene symbols on ALL 4 models
# ═════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("FIXING GENE SYMBOLS ON ALL MODELS")
print("=" * 65)

mapping['EntrezID'] = mapping['EntrezID'].astype(str)

for model_num in range(1, 5):
    result_path = BASE / 'analysis_v2' / f'results_model{model_num}' / 'DE_results.csv'
    if not result_path.exists():
        print(f"\n  WARNING: {result_path} not found, skipping")
        continue

    res = pd.read_csv(result_path, index_col=0)

    # Index is Entrez IDs — clean any .0 suffix
    res.index = res.index.astype(str).str.replace(r'\.0$', '', regex=True)
    res.index.name = 'EntrezID'

    # Merge with symbol mapping
    res_reset = res.reset_index()
    res_merged = res_reset.merge(mapping, on='EntrezID', how='left')

    # Use symbol where available, fall back to Entrez ID
    res_merged['GeneSymbol'] = res_merged['Symbol'].fillna(res_merged['EntrezID'])
    res_merged = res_merged.set_index('GeneSymbol')
    res_merged = res_merged.drop(columns=['EntrezID', 'Symbol'], errors='ignore')

    # Save (overwrite with symbol-indexed version)
    res_merged.to_csv(result_path)

    n_mapped = res_merged.index.isin(mapping['Symbol'].values).sum()
    print(f"\n  Model {model_num}: {n_mapped:,}/{len(res_merged):,} genes mapped to symbols")

# ═════════════════════════════════════════════════════════════
# Final summary
# ═════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("MODEL PROGRESSION SUMMARY")
print("=" * 65)

# Re-read all models to get consistent counts
for model_num in range(1, 5):
    res = pd.read_csv(BASE / 'analysis_v2' / f'results_model{model_num}' / 'DE_results.csv', index_col=0)
    n_sig = (res['padj'] < 0.05).sum()
    pct = n_sig / len(res) * 100
    labels = {
        1: "~ Diagnosis",
        2: "~ Age + Sex + Diagnosis",
        3: "~ Age + Sex + Braak + Diagnosis",
        4: "~ Age + Sex + Braak + CellScores + Diagnosis",
    }
    print(f"  Model {model_num}: {n_sig:>6,} sig genes ({pct:5.1f}%)  {labels[model_num]}")

# Quick check on Model 4 top genes
res4 = pd.read_csv(BASE / 'analysis_v2' / 'results_model4' / 'DE_results.csv', index_col=0)
sig4 = res4[res4['padj'] < 0.05].sort_values('padj')

print(f"\nModel 4 — Top 10 upregulated:")
for gene, row in sig4[sig4['log2FoldChange'] > 0].head(10).iterrows():
    print(f"  {str(gene):<15} LFC={row['log2FoldChange']:7.3f}  padj={row['padj']:.2e}")

print(f"\nModel 4 — Top 10 downregulated:")
for gene, row in sig4[sig4['log2FoldChange'] < 0].head(10).iterrows():
    print(f"  {str(gene):<15} LFC={row['log2FoldChange']:7.3f}  padj={row['padj']:.2e}")

# Check key genes
for gene in ['APOE', 'TREM2', 'ABCA7', 'EGR1', 'EGR2', 'FOSB', 'CCL4']:
    if gene in res4.index:
        row = res4.loc[gene]
        sig_str = "SIG" if row['padj'] < 0.05 else "ns"
        print(f"  {gene:<10} LFC={row['log2FoldChange']:7.3f}  padj={row['padj']:.2e}  [{sig_str}]")

print("\n" + "=" * 65)
print("DONE.")
print("  All 4 model results saved with gene symbols.")
print("  Next: python scripts/06_enrichment_and_verification.py")
print("=" * 65)
