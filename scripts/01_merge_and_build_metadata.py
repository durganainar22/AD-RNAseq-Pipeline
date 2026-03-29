#!/usr/bin/env python
"""
01_merge_and_build_metadata.py
-------------------------------
Merge 289 individual GSM count files into one count matrix,
parse the GEO SOFT file for clinical metadata, and produce
a single clean metadata table.

Inputs:
    data/raw/GSM*.tsv.gz                 (289 per-sample count files)
    data/raw/GSE125583_family.soft.gz    (GEO SOFT annotation)

Outputs:
    data/raw/merged_counts.tsv           (20,024 genes x 289 samples)
    data/metadata/sample_metadata_full.csv
        Columns: SampleID, GEO_ID, Diagnosis, Age, Sex, Sex_binary,
                 BraakScore, BraakScore_numeric
"""

import gzip
import glob
import os
import sys
import pandas as pd
from pathlib import Path

BASE = Path('path here')
os.chdir(BASE)

print("=" * 65)
print("SCRIPT 1: Merge Counts + Build Metadata")
print("Dataset: GSE125583 — Fusiform Gyrus Bulk RNA-seq")
print("=" * 65)

# ═════════════════════════════════════════════════════════════
# PART A: Merge 289 sample files into count matrix
# ═════════════════════════════════════════════════════════════
print("\n--- PART A: Merging sample count files ---\n")

sample_files = sorted(glob.glob('data/raw/GSM*.tsv.gz'))

if len(sample_files) == 0:
    print("ERROR: No GSM*.tsv.gz files found in data/raw/")
    sys.exit(1)

print(f"Found {len(sample_files)} sample files")

# Read first file to get gene IDs and inspect structure
first = pd.read_csv(sample_files[0], sep='\t', comment='#')
print(f"Columns per file: {first.columns.tolist()}")
print(f"Genes per file:   {len(first):,}")

# Build count matrix: genes (rows) x samples (columns)
count_matrix = pd.DataFrame(index=first['ID_REF'])
count_matrix.index.name = 'Gene'

for i, filepath in enumerate(sample_files):
    sample_name = os.path.basename(filepath).replace('.tsv.gz', '')
    df = pd.read_csv(filepath, sep='\t', comment='#')
    count_matrix[sample_name] = df['count'].values

    if (i + 1) % 50 == 0 or (i + 1) == len(sample_files):
        print(f"  Merged {i + 1}/{len(sample_files)} samples")

print(f"\nCount matrix: {count_matrix.shape[0]:,} genes x {count_matrix.shape[1]} samples")

# Sanity checks
lib_sizes = count_matrix.sum(axis=0)
zero_genes = (count_matrix == 0).all(axis=1).sum()
print(f"Library size range: {lib_sizes.min():,.0f} – {lib_sizes.max():,.0f}")
print(f"Genes with all zeros: {zero_genes}")

# Save
out_counts = BASE / 'data' / 'raw' / 'merged_counts.tsv'
count_matrix.to_csv(out_counts, sep='\t')
print(f"Saved: {out_counts}")

# ═════════════════════════════════════════════════════════════
# PART B: Parse SOFT file for clinical metadata
# ═════════════════════════════════════════════════════════════
print("\n--- PART B: Extracting metadata from SOFT file ---\n")

soft_path = BASE / 'data' / 'raw' / 'GSE125583_family.soft.gz'

if not soft_path.exists():
    print(f"ERROR: {soft_path} not found.")
    sys.exit(1)

records = {}
gsm = None
chars = {}

with gzip.open(soft_path, 'rt') as f:
    for line in f:
        if line.startswith('^SAMPLE'):
            if gsm:
                records[gsm] = chars.copy()
            gsm = line.strip().split(' = ')[1]
            chars = {}
        elif line.startswith('!Sample_characteristics_ch1'):
            val = line.strip().split(' = ', 1)[1]
            if ':' in val:
                k, v = val.split(':', 1)
                chars[k.strip()] = v.strip()
    if gsm:
        records[gsm] = chars

soft_df = pd.DataFrame(records).T
soft_df.index.name = 'GEO_ID'
soft_df = soft_df.reset_index()

print(f"Parsed {len(soft_df)} samples from SOFT file")
print(f"Fields found: {soft_df.columns.tolist()}")

# ═════════════════════════════════════════════════════════════
# PART C: Build clean metadata table
# ═════════════════════════════════════════════════════════════
print("\n--- PART C: Building metadata table ---\n")

# Create base metadata from count matrix column names
meta = pd.DataFrame({
    'SampleID': count_matrix.columns,
    'GEO_ID': [name.split('_')[0] for name in count_matrix.columns],
})

# Merge with SOFT data
merged = meta.merge(soft_df, on='GEO_ID', how='left')

# Standardize column names
rename_map = {
    'age': 'Age',
    'Sex': 'Sex',
    'braak.score': 'BraakScore',
    'diagnosis': 'Diagnosis_GEO',
}
for old, new in rename_map.items():
    if old in merged.columns:
        merged = merged.rename(columns={old: new})

# --- Diagnosis ---
if 'Diagnosis_GEO' in merged.columns:
    diag_map = {
        "Alzheimer's disease": 'AD',
        "Alzheimer's Disease": 'AD',
        'AD': 'AD',
        'Control': 'Control',
        'control': 'Control',
    }
    merged['Diagnosis'] = merged['Diagnosis_GEO'].map(diag_map)
    unmapped = merged['Diagnosis'].isna().sum()
    if unmapped > 0:
        print(f"WARNING: {unmapped} samples with unmapped diagnosis:")
        print(merged[merged['Diagnosis'].isna()]['Diagnosis_GEO'].unique())

# --- Age (numeric) ---
if 'Age' in merged.columns:
    merged['Age'] = pd.to_numeric(merged['Age'], errors='coerce')

# --- Sex (binary: M=1, F=0) ---
if 'Sex' in merged.columns:
    merged['Sex_binary'] = merged['Sex'].map({'M': 1, 'F': 0})

# --- Braak score (Roman numerals -> integer) ---
braak_map = {'0': 0, 'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6}
if 'BraakScore' in merged.columns:
    merged['BraakScore_numeric'] = (
        merged['BraakScore'].astype(str).str.strip().map(braak_map)
    )

# ═════════════════════════════════════════════════════════════
# PART D: Report and save
# ═════════════════════════════════════════════════════════════
print("=" * 45)
print("METADATA SUMMARY")
print("=" * 45)

print(f"\nTotal samples: {len(merged)}")
print(f"\nDiagnosis:")
print(merged['Diagnosis'].value_counts().to_string())

if 'Age' in merged.columns:
    print(f"\nAge by diagnosis:")
    print(merged.groupby('Diagnosis')['Age'].agg(
        ['count', 'mean', 'std', 'min', 'max']).round(1).to_string())

if 'Sex' in merged.columns:
    print(f"\nSex by diagnosis:")
    print(pd.crosstab(merged['Diagnosis'], merged['Sex']).to_string())

if 'BraakScore' in merged.columns:
    print(f"\nBraak distribution:")
    print(merged['BraakScore'].value_counts().sort_index().to_string())

# Missing values
print(f"\nMissing values:")
check_cols = [c for c in ['Age', 'Sex', 'Sex_binary', 'BraakScore',
                           'BraakScore_numeric', 'Diagnosis']
              if c in merged.columns]
print(merged[check_cols].isna().sum().to_string())

# Save
(BASE / 'data' / 'metadata').mkdir(parents=True, exist_ok=True)
out_meta = BASE / 'data' / 'metadata' / 'sample_metadata_full.csv'
merged.to_csv(out_meta, index=False)
print(f"\nSaved: {out_meta}")
print(f"Columns: {merged.columns.tolist()}")

print("\n" + "=" * 65)
print("DONE.")
print(f"  Count matrix:  data/raw/merged_counts.tsv ({count_matrix.shape[0]:,} x {count_matrix.shape[1]})")
print(f"  Metadata:      data/metadata/sample_metadata_full.csv ({len(merged)} samples)")
print("  Next: python scripts/02_qc_and_filter.py")
print("=" * 65)
