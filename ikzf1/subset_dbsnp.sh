#!/usr/bin/bash

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO hg19
DB19=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz 
OUT19=./data/dbsnp_ikzf1_hg19.vcf

# CMD
bcftools view -r NC_000007.13:49300000-51310000 $DB19 |\
    bcftools annotate -x INFO -O v -o $OUT19 -

# Define IO hg38
DB38=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39.gz
OUT38=./data/dbsnp_ikzf1_hg38.vcf

# CMD
bcftools view -r NC_000007.14:49200000-51300000 $DB38 |\
    bcftools annotate -x INFO -O v -o $OUT38 -

