#!/bin/bash

###########################################
# 16S Amplicon Data Re-analysis Pipeline  #
#                                         #
# This script performs re-analysis of     #
# public 16S rRNA sequencing data from SRA.
# It covers the full workflow from:       #
#   1. Downloading FASTQ files via SRR    #
#   2. Generating QIIME2-compatible input #
#   3. Denoising using DADA2              #
#   4. Taxonomic classification           #
#   5. Exporting OTU tables               #
#                                         #
# To use this script:                     #
#   - Provide a file named `srr.txt`      #
#     with one SRR accession per line     #
#   - Modify the user configuration block #
#     below to match your environment     #
###########################################

##############################
#        User Settings       #
##############################

# Tool paths
SRATOOLKIT_PATH=/home/han_xiao/xiaohan/Software/sratoolkit.3.2.0-centos_linux64/bin
QIIME_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/qiime2-amplicon-2024.10/bin/qiime
BIOM_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/qiime2-amplicon-2024.10/bin/biom
RUSH_PATH=/home/han_xiao/xiaohan/Software/miniconda3/envs/rush/bin/rush

# Working directory for all outputs
WORK_DIR=/home/han_xiao/xiaohan/Others/Ruwen/SRP542009/01.data
SRR_LIST=$WORK_DIR/srr.txt   # List of SRR IDs, one per line

# Path to pre-trained classifier
CLASSIFIER_DB=/home/han_xiao/xiaohan/Database/16S/gg-13-8-99-nb-classifier.qza

CPU=40

##############################
#         Main Steps         #
##############################

# Step 1: Download and convert SRR data
echo Step 1 downloading and converting SRR data to FASTQ
mkdir -p "$WORK_DIR" && cd "$WORK_DIR"

if [ ! -f "$SRR_LIST" ]; then
    echo Error SRR list file not found
    exit 1
fi

cat "$SRR_LIST" | "$RUSH_PATH" -j 40 -k '
    '"$SRATOOLKIT_PATH"'/prefetch {}
    '"$SRATOOLKIT_PATH"'/fasterq-dump --split-3 --outdir ./{}/ ./{}/{}.sra
    rm -f ./{}/{}.sra
'

# Step 2: Generate QIIME2 manifest file
echo Step 2 generating QIIME2 manifest file
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

# Step 3: Import FASTQ files to QIIME2 .qza format
echo Step 3 importing FASTQ into QIIME2 format
"$QIIME_PATH" tools import \
  --input-format PairedEndFastqManifestPhred33V2 \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path data.qza

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

CLASSIFIER_DB=/home/han_xiao/xiaohan/Database/16S/gg-13-8-99-nb-classifier.qza

"$QIIME_PATH" feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_DB" \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs $CPU

mkdir -p otu/table otu/taxonomy

"$QIIME_PATH" tools export --input-path table.qza --output-path otu/table
"$QIIME_PATH" tools export --input-path taxonomy.qza --output-path otu/taxonomy

"$BIOM_PATH" convert -i otu/table/feature-table.biom -o otu/otu_table.tsv --to-tsv
echo Pipeline completed Results are saved in $WORK_DIR/otu

