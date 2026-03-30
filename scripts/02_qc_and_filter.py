#!/usr/bin/env python
"""
02_qc_and_filter.py
--------------------
Quality control checks and low-count gene filtering.

Inputs:
    data/raw/merged_counts.tsv                (30,727 genes x 289 samples)
    data/metadata/sample_metadata_full.csv    (289 samples with clinical info)

Outputs:
    data/processed/02_filtered_counts.csv     (filtered genes x 289 samples)
    results/qc/01_library_sizes.pdf
    results/qc/02_gene_detection.pdf
    results/qc/03_count_distribution.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

BASE = Path('/path to project')
os.chdir(BASE)

print("=" * 65)
print("SCRIPT 2: Quality Control and Filtering")
print("=" * 65)

# ═════════════════════════════════════════════════════════════
# Load data
# ═════════════════════════════════════════════════════════════
print("\n--- Loading data ---\n")

counts_file = BASE / 'data' / 'raw' / 'merged_counts.tsv'
meta_file = BASE / 'data' / 'metadata' / 'sample_metadata_full.csv'

for f in [counts_file, meta_file]:
    if not f.exists():
        print(f"ERROR: {f} not found. Run 01_merge_and_build_metadata.py first.")
        sys.exit(1)

counts = pd.read_csv(counts_file, sep='\t', index_col=0)
meta = pd.read_csv(meta_file)

print(f"Count matrix: {counts.shape[0]:,} genes x {counts.shape[1]} samples")
print(f"Metadata:     {len(meta)} samples")

# ═════════════════════════════════════════════════════════════
# PART A: Library size QC
# ═════════════════════════════════════════════════════════════
print("\n--- PART A: Library size QC ---\n")

(BASE / 'results' / 'qc').mkdir(parents=True, exist_ok=True)

lib_sizes = counts.sum(axis=0) / 1e6  # millions of reads

print(f"Library size (millions of reads):")
print(f"  Min:    {lib_sizes.min():.2f}M")
print(f"  Median: {lib_sizes.median():.2f}M")
print(f"  Max:    {lib_sizes.max():.2f}M")
print(f"  Mean:   {lib_sizes.mean():.2f}M")

# Flag potential outliers (< 1M reads)
low_lib = lib_sizes[lib_sizes < 1.0]
if len(low_lib) > 0:
    print(f"\n  WARNING: {len(low_lib)} samples with < 1M reads:")
    for name, size in low_lib.items():
        print(f"    {name}: {size:.2f}M")
else:
    print(f"  No low-library-size outliers detected.")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Histogram
axes[0].hist(lib_sizes, bins=40, color='steelblue', edgecolor='white')
axes[0].axvline(lib_sizes.median(), color='red', ls='--', lw=1.5,
                label=f'Median: {lib_sizes.median():.1f}M')
axes[0].set_xlabel('Library Size (Millions of Reads)')
axes[0].set_ylabel('Number of Samples')
axes[0].set_title('Library Size Distribution')
axes[0].legend()

# Sorted bar plot
lib_sorted = lib_sizes.sort_values()
colors = ['steelblue'] * len(lib_sorted)
axes[1].bar(range(len(lib_sorted)), lib_sorted.values, color=colors, width=1.0)
axes[1].axhline(lib_sizes.median(), color='red', ls='--', lw=1.5)
axes[1].set_xlabel('Samples (sorted)')
axes[1].set_ylabel('Library Size (Millions)')
axes[1].set_title('Library Sizes (sorted)')

plt.tight_layout()
plt.savefig(BASE / 'results' / 'qc' / '01_library_sizes.pdf', dpi=300)
plt.close()
print("  Saved: results/qc/01_library_sizes.pdf")

# ═════════════════════════════════════════════════════════════
# PART B: Gene detection QC
# ═════════════════════════════════════════════════════════════
print("\n--- PART B: Gene detection QC ---\n")

genes_detected = (counts > 0).sum(axis=0)

print(f"Genes detected per sample:")
print(f"  Min:    {genes_detected.min():,}")
print(f"  Median: {int(genes_detected.median()):,}")
print(f"  Max:    {genes_detected.max():,}")

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(genes_detected, bins=40, color='steelblue', edgecolor='white')
ax.axvline(genes_detected.median(), color='red', ls='--', lw=1.5,
           label=f'Median: {int(genes_detected.median()):,}')
ax.set_xlabel('Number of Genes Detected (count > 0)')
ax.set_ylabel('Number of Samples')
ax.set_title('Gene Detection per Sample')
ax.legend()
plt.tight_layout()
plt.savefig(BASE / 'results' / 'qc' / '02_gene_detection.pdf', dpi=300)
plt.close()
print("  Saved: results/qc/02_gene_detection.pdf")

# ═════════════════════════════════════════════════════════════
# PART C: Count distribution (before filtering)
# ═════════════════════════════════════════════════════════════
print("\n--- PART C: Count distribution ---\n")

# Log2 of mean counts per gene (across all samples)
mean_counts = counts.mean(axis=1)
log2_mean = np.log2(mean_counts + 1)

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(log2_mean, bins=80, color='steelblue', edgecolor='white')
ax.axvline(np.log2(10 + 1), color='red', ls='--', lw=1.5,
           label='count = 10 (filter threshold)')
ax.set_xlabel('Log2(Mean Count + 1)')
ax.set_ylabel('Number of Genes')
ax.set_title('Gene Expression Distribution (before filtering)')
ax.legend()
plt.tight_layout()
plt.savefig(BASE / 'results' / 'qc' / '03_count_distribution.pdf', dpi=300)
plt.close()
print("  Saved: results/qc/03_count_distribution.pdf")

# ═════════════════════════════════════════════════════════════
# PART D: Filter low-count genes
# ═════════════════════════════════════════════════════════════
print("\n--- PART D: Filtering low-count genes ---\n")

min_count = 10
min_fraction = 0.10
min_samples = int(min_fraction * counts.shape[1])

print(f"Filter criteria:")
print(f"  Keep genes with count >= {min_count} in >= {min_samples} samples")
print(f"  ({min_fraction*100:.0f}% of {counts.shape[1]} samples)")

genes_pass = (counts >= min_count).sum(axis=1) >= min_samples
counts_filtered = counts[genes_pass]

n_removed = counts.shape[0] - counts_filtered.shape[0]
pct_kept = counts_filtered.shape[0] / counts.shape[0] * 100

print(f"\nFiltering results:")
print(f"  Before:   {counts.shape[0]:,} genes")
print(f"  After:    {counts_filtered.shape[0]:,} genes")
print(f"  Removed:  {n_removed:,} genes ({100 - pct_kept:.1f}%)")
print(f"  Retained: {pct_kept:.1f}%")

# ═════════════════════════════════════════════════════════════
# PART E: Verify sample-metadata alignment and save
# ═════════════════════════════════════════════════════════════
print("\n--- PART E: Final checks and save ---\n")

# Match samples between counts and metadata
count_samples = set(counts_filtered.columns)
meta_samples = set(meta['SampleID'])

common = sorted(count_samples & meta_samples)
only_counts = count_samples - meta_samples
only_meta = meta_samples - count_samples

print(f"Sample alignment:")
print(f"  In counts:   {len(count_samples)}")
print(f"  In metadata: {len(meta_samples)}")
print(f"  Common:      {len(common)}")
if only_counts:
    print(f"  Only in counts (no metadata): {len(only_counts)}")
if only_meta:
    print(f"  Only in metadata (no counts): {len(only_meta)}")

# Keep only common samples, aligned
counts_filtered = counts_filtered[common]
meta_filtered = meta[meta['SampleID'].isin(common)].copy()
meta_filtered = meta_filtered.set_index('SampleID').loc[common].reset_index()

print(f"\nFinal dimensions:")
print(f"  Genes:    {counts_filtered.shape[0]:,}")
print(f"  Samples:  {counts_filtered.shape[1]}")
print(f"  AD:       {(meta_filtered['Diagnosis'] == 'AD').sum()}")
print(f"  Control:  {(meta_filtered['Diagnosis'] == 'Control').sum()}")

# Save
(BASE / 'data' / 'processed').mkdir(parents=True, exist_ok=True)

out_counts = BASE / 'data' / 'processed' / '02_filtered_counts.csv'
counts_filtered.to_csv(out_counts)
print(f"\nSaved: {out_counts}")

out_meta = BASE / 'data' / 'processed' / '02_filtered_metadata.csv'
meta_filtered.to_csv(out_meta, index=False)
print(f"Saved: {out_meta}")

print("\n" + "=" * 65)
print("DONE.")
print(f"  Filtered counts: {counts_filtered.shape[0]:,} genes x {counts_filtered.shape[1]} samples")
print(f"  QC plots:        results/qc/")
print("  Next: python scripts/03_build_covariates.py")
print("=" * 65)
