#!/bin/bash

echo "=== Downloading GSE125583 from GEO ==="
echo "Dataset: AD Fusiform Gyrus (219 AD + 70 Controls)"
echo ""

cd *project directory path here*

# Download the RAW tar file (contains all sample data)
echo "Downloading RAW data (87.4 MB)..."
wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE125583&format=file" -O GSE125583_RAW.tar

echo "✓ Download complete!"
echo ""

# Extract the tar file
echo "Extracting files..."
tar -xf GSE125583_RAW.tar

echo "✓ Extraction complete!"
echo ""

# List what we got
echo "Downloaded files:"
ls -lh *.tsv.gz 2>/dev/null | head -10

echo ""
echo "Total samples:"
ls -1 *.tsv.gz 2>/dev/null | wc -l

echo ""
echo "✓ Data ready for analysis!"
