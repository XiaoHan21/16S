#!/bin/bash

# ===============================================================
# 16S Amplicon Data Re-analysis Pipeline
#
# This script performs re-analysis of public 16S rRNA sequencing 
# data from SRA using QIIME2.
#
# Steps:
#   1. Download FASTQ files from SRR list
#   2. Generate QIIME2-compatible manifest file
#   3. Import data to QIIME2
#   4. Denoise with DADA2
#   5. Taxonomic classification with sklearn classifier
#   6. Export OTU tables
#
# Requirements:
#   - srr.txt: file with one SRR ID per line
#   - Set correct paths to software and database below
# ===============================================================

# ---------------- User Configuration ---------------- #

SRATOOLKIT_PATH=/home/han_xiao/xiaohan/Software/sratoolkit.3.2.0-centos_linux64/bin
QIIME_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/qiime2-amplicon-2024.10/bin/qiime
BIOM_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/qiime2-amplicon-2024.10/bin/biom
RUSH_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/rush/bin/rush

WORK_DIR=/home/han_xiao/xiaohan/Others/Ruwen/SRP542009/01.data
SRR_LIST=$WORK_DIR/srr.txt
CLASSIFIER_DB=/home/han_xiao/xiaohan/Database/16S/gg-13-8-99-nb-classifier.qza
CPU=40

# ---------------- Step 1: Download SRR FASTQ ---------------- #

echo Step 1: downloading and converting SRR data to FASTQ
mkdir -p "$WORK_DIR" && cd "$WORK_DIR"

if [ ! -f "$SRR_LIST" ]; then
    echo Error: SRR list file not found: $SRR_LIST
    exit 1
fi

cat "$SRR_LIST" | "$RUSH_PATH" -j $CPU -k '
    '"$SRATOOLKIT_PATH"'/prefetch {}
    '"$SRATOOLKIT_PATH"'/fasterq-dump --split-3 --outdir ./{}/ ./{}/{}.sra
    rm -f ./{}/{}.sra
'

# ---------------- Step 2: Generate manifest ---------------- #

echo Step 2: generating QIIME2 manifest file
MANIFEST_FILE="$WORK_DIR/samples.manifest"

{
  echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"
  base_dir=$(pwd)
  for sample_dir in SRR*; do
      if [ -d "$sample_dir" ]; then
          sample_id="$sample_dir"
          fwd_file="$base_dir/$sample_dir/${sample_id}_1.fastq"
          rev_file="$base_dir/$sample_dir/${sample_id}_2.fastq"
          echo -e "${sample_id}\t${fwd_file}\t${rev_file}"
      fi
  done
} > "$MANIFEST_FILE"

# ---------------- Step 3: Import to QIIME2 ---------------- #

echo Step 3: importing FASTQ files to QIIME2
"$QIIME_PATH" tools import \
  --input-format PairedEndFastqManifestPhred33V2 \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path data.qza

# ---------------- Step 4: DADA2 denoising ---------------- #

echo Step 4: denoising using DADA2
"$QIIME_PATH" dada2 denoise-paired \
  --i-demultiplexed-seqs data.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza \
  --p-n-threads $CPU

# ---------------- Step 5: Taxonomic classification ---------------- #

echo Step 5: classifying sequences
"$QIIME_PATH" feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_DB" \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs $CPU

# ---------------- Step 6: Export results ---------------- #

echo Step 6: exporting OTU results
mkdir -p otu/table otu/taxonomy

"$QIIME_PATH" tools export --input-path table.qza --output-path otu/table
"$QIIME_PATH" tools export --input-path taxonomy.qza --output-path otu/taxonomy

"$BIOM_PATH" convert -i otu/table/feature-table.biom -o otu/otu_table.tsv --to-tsv

echo Pipeline completed. Results saved in: $WORK_DIR/otu
