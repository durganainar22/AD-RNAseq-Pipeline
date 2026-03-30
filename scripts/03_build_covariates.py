#!/usr/bin/env python
"""
03_build_covariates.py
-----------------------
Build Entrez→Symbol mapping, normalize counts, convert to gene symbols,
compute cell-type deconvolution scores, and center continuous covariates.

Inputs:
    data/processed/02_filtered_counts.csv          (20,024 genes x 289 samples)
    data/metadata/sample_metadata_full.csv         (clinical metadata)
    data/annotations/Homo_sapiens.gene_info.gz     (NCBI gene annotation)

Outputs:
    data/annotations/entrez_to_symbol.csv          (Entrez ID → gene symbol map)
    data/processed/04_deseq2_normalized_counts.csv (normalized, Entrez IDs)
    data/processed/04_deseq2_normalized_counts_symbols.csv (normalized, gene symbols)
    data/metadata/sample_metadata_covariates.csv   (final metadata for DESeq2)
"""

import gzip
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

BASE = Path('/path to project')
os.chdir(BASE)

print("=" * 65)
print("SCRIPT 3: Build Covariates")
print("  - Entrez → Symbol mapping")
print("  - Library-size normalization")
print("  - Cell-type deconvolution scores")
print("  - Centered covariates")
print("=" * 65)

# ═════════════════════════════════════════════════════════════
# Load data
# ═════════════════════════════════════════════════════════════
print("\n--- Loading data ---\n")

counts = pd.read_csv(BASE / 'data' / 'processed' / '02_filtered_counts.csv', index_col=0)
meta = pd.read_csv(BASE / 'data' / 'metadata' / 'sample_metadata_full.csv')

gene_info_path = BASE / 'data' / 'annotations' / 'Homo_sapiens.gene_info.gz'
if not gene_info_path.exists():
    print(f"ERROR: {gene_info_path} not found.")
    print("Download: wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
    print("Place in: data/annotations/")
    sys.exit(1)

print(f"Filtered counts: {counts.shape[0]:,} genes x {counts.shape[1]} samples")
print(f"Metadata:        {len(meta)} samples")

# ═════════════════════════════════════════════════════════════
# PART A: Build Entrez → Symbol mapping
# ═════════════════════════════════════════════════════════════
print("\n--- PART A: Building Entrez → Symbol mapping ---\n")

(BASE / 'data' / 'annotations').mkdir(parents=True, exist_ok=True)

rows = []
with gzip.open(gene_info_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        entrez = parts[1]
        symbol = parts[2]
        rows.append((entrez, symbol))

mapping = pd.DataFrame(rows, columns=['EntrezID', 'Symbol']).drop_duplicates()

out_map = BASE / 'data' / 'annotations' / 'entrez_to_symbol.csv'
mapping.to_csv(out_map, index=False)
print(f"Total mappings: {len(mapping):,}")
print(f"Saved: {out_map}")

# Check how many of our genes have a mapping
counts.index = counts.index.astype(str)
mapping['EntrezID'] = mapping['EntrezID'].astype(str)

mapped_genes = set(counts.index) & set(mapping['EntrezID'])
unmapped_genes = set(counts.index) - set(mapping['EntrezID'])
print(f"Genes in count matrix: {counts.shape[0]:,}")
print(f"  With symbol:  {len(mapped_genes):,}")
print(f"  No symbol:    {len(unmapped_genes):,}")

# ═════════════════════════════════════════════════════════════
# PART B: Library-size normalization
# ═════════════════════════════════════════════════════════════
print("\n--- PART B: Library-size normalization ---\n")

# Simple median-ratio normalization (similar to DESeq2 size factors)
# For cell-type scoring we just need CPM-style normalization
lib_sizes = counts.sum(axis=0)
size_factors = lib_sizes / lib_sizes.median()
normalized = counts.div(size_factors, axis=1)

out_norm = BASE / 'data' / 'processed' / '04_deseq2_normalized_counts.csv'
normalized.to_csv(out_norm)
print(f"Normalized counts saved: {out_norm}")
print(f"  Shape: {normalized.shape[0]:,} x {normalized.shape[1]}")

# ═════════════════════════════════════════════════════════════
# PART C: Convert normalized counts to gene symbols
# ═════════════════════════════════════════════════════════════
print("\n--- PART C: Converting to gene symbols ---\n")

# Merge normalized counts with mapping
norm_with_symbol = normalized.copy()
norm_with_symbol.index.name = 'EntrezID'
norm_with_symbol = norm_with_symbol.reset_index()
norm_with_symbol['EntrezID'] = norm_with_symbol['EntrezID'].astype(str)

merged = norm_with_symbol.merge(mapping, on='EntrezID', how='left')

# Drop genes without symbol, group duplicates by mean
mapped = merged.dropna(subset=['Symbol']).copy()
mapped = mapped.drop(columns=['EntrezID'])
mapped = mapped.groupby('Symbol').mean()

out_symbols = BASE / 'data' / 'processed' / '04_deseq2_normalized_counts_symbols.csv'
mapped.to_csv(out_symbols)
print(f"Symbol-named matrix: {mapped.shape[0]:,} genes x {mapped.shape[1]} samples")
print(f"Saved: {out_symbols}")

# Quick check for key genes
key_genes = ['APOE', 'TREM2', 'GFAP', 'RBFOX3', 'MBP', 'AIF1']
found = [g for g in key_genes if g in mapped.index]
missing = [g for g in key_genes if g not in mapped.index]
print(f"Key gene check: {len(found)}/{len(key_genes)} found")
if missing:
    print(f"  Missing: {missing}")

# ═════════════════════════════════════════════════════════════
# PART D: Compute cell-type deconvolution scores
# ═════════════════════════════════════════════════════════════
print("\n--- PART D: Cell-type deconvolution scores ---\n")

# Log2-transform for scoring (add 1 to avoid log(0))
log_norm = np.log2(mapped + 1)

# Marker gene lists (well-established brain cell-type markers)
MARKERS = {
    'MicrogliaScore': ['TYROBP', 'AIF1', 'CSF1R', 'SPI1', 'TREM2',
                       'C1QA', 'C1QB', 'C1QC'],
    'NeuronScore':    ['RBFOX3', 'SYT1', 'SLC17A7', 'SNAP25', 'MAP2',
                       'NEFL', 'NEFM'],
    'AstrocyteScore': ['GFAP', 'ALDH1L1', 'AQP4', 'SLC1A3', 'S100B'],
    'OligoScore':     ['MBP', 'MOG', 'PLP1', 'MAG', 'OLIG1', 'OLIG2'],
}

# Align metadata to count matrix sample order
meta = meta.set_index('SampleID')
common = sorted(set(meta.index) & set(mapped.columns))
meta = meta.loc[common].copy()
log_norm = log_norm[common]

for score_name, gene_list in MARKERS.items():
    present = [g for g in gene_list if g in log_norm.index]
    if len(present) == 0:
        print(f"  WARNING: No genes found for {score_name}")
        meta[score_name] = np.nan
    else:
        meta[score_name] = log_norm.loc[present].mean(axis=0).values
        print(f"  {score_name}: {len(present)}/{len(gene_list)} markers used: {present}")

# ═════════════════════════════════════════════════════════════
# PART E: Center continuous covariates
# ═════════════════════════════════════════════════════════════
print("\n--- PART E: Centering covariates ---\n")

# Ensure numeric
meta['Age'] = pd.to_numeric(meta['Age'], errors='coerce')
meta['Sex_binary'] = pd.to_numeric(meta['Sex_binary'], errors='coerce')
meta['BraakScore_numeric'] = pd.to_numeric(meta['BraakScore_numeric'], errors='coerce')

# Center
meta['Age_c'] = meta['Age'] - meta['Age'].mean()
meta['Braak_c'] = meta['BraakScore_numeric'] - meta['BraakScore_numeric'].mean()

print(f"Age:   mean={meta['Age'].mean():.1f}, centered range: "
      f"[{meta['Age_c'].min():.1f}, {meta['Age_c'].max():.1f}]")
print(f"Braak: mean={meta['BraakScore_numeric'].mean():.2f}, centered range: "
      f"[{meta['Braak_c'].min():.2f}, {meta['Braak_c'].max():.2f}]")

# ═════════════════════════════════════════════════════════════
# PART F: Report and save
# ═════════════════════════════════════════════════════════════
print("\n--- PART F: Final report and save ---\n")

# Missing values check
print("Missing values in covariate columns:")
cov_cols = ['Diagnosis', 'Age', 'Age_c', 'Sex_binary', 'BraakScore_numeric',
            'Braak_c', 'MicrogliaScore', 'NeuronScore', 'AstrocyteScore', 'OligoScore']
for col in cov_cols:
    if col in meta.columns:
        n_miss = meta[col].isna().sum()
        if n_miss > 0:
            print(f"  {col}: {n_miss} missing")

print(f"\nTotal samples: {len(meta)}")
print(f"  AD:      {(meta['Diagnosis'] == 'AD').sum()}")
print(f"  Control: {(meta['Diagnosis'] == 'Control').sum()}")

# Cell score summary by diagnosis
print(f"\nCell-type scores by diagnosis:")
score_cols = ['MicrogliaScore', 'NeuronScore', 'AstrocyteScore', 'OligoScore']
for col in score_cols:
    if col in meta.columns and not meta[col].isna().all():
        ad_mean = meta[meta['Diagnosis'] == 'AD'][col].mean()
        ctl_mean = meta[meta['Diagnosis'] == 'Control'][col].mean()
        direction = "AD > Ctl" if ad_mean > ctl_mean else "Ctl > AD"
        print(f"  {col:<18} AD={ad_mean:.2f}  Ctl={ctl_mean:.2f}  ({direction})")

# Save
meta = meta.reset_index()
out_cov = BASE / 'data' / 'metadata' / 'sample_metadata_covariates.csv'
meta.to_csv(out_cov, index=False)
print(f"\nSaved: {out_cov}")

# Preview
print(f"\nPreview of covariate columns:")
preview_cols = ['SampleID', 'Diagnosis', 'Age_c', 'Sex_binary', 'Braak_c',
                'MicrogliaScore', 'NeuronScore', 'AstrocyteScore', 'OligoScore']
preview_cols = [c for c in preview_cols if c in meta.columns]
print(meta[preview_cols].head(5).to_string())

print("\n" + "=" * 65)
print("DONE.")
print(f"  Entrez→Symbol map:  data/annotations/entrez_to_symbol.csv ({len(mapping):,} mappings)")
print(f"  Normalized counts:  data/processed/04_deseq2_normalized_counts_symbols.csv")
print(f"  Covariates:         data/metadata/sample_metadata_covariates.csv ({len(meta)} samples)")
print("  Next: python scripts/04_baseline_deseq2.py")
print("=" * 65)
