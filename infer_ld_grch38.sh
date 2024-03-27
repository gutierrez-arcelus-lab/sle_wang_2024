#!/usr/bin/sh

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

IDS=./data/eur_kgp.txt
REGION=$( awk '{ print $1 }' ./data/ikzf1_coords_grch38.txt)
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/GRCh38/Genotype_VCFs/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.vcf.gz
PREFIX=data/chr${REGION}_grch38
VCFOUT=${PREFIX}.vcf.gz

bcftools view --threads 4 -r chr${REGION} $VCF1K |\
    bcftools view --threads 4 --samples-file $IDS --force-samples |\
    bcftools norm --threads 4 -m - - |\
    bcftools view --threads 4 --min-ac 2:minor - |\
    bcftools view --threads 4 -e INFO/AC==INFO/AN - |\
    bcftools annotate --threads 4 -x FORMAT -O z -o $VCFOUT -

plink --vcf $VCFOUT --r2 square spaces --real-ref-alleles --threads 4 --out ${PREFIX}
