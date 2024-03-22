#!/usr/bin/sh

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

IDS=ikzf1/data/eur_kgp.txt
REGION=1:205643772-207643772
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
PREFIX=data/chr${REGION}
VCFOUT=${PREFIX}.vcf.gz

bcftools view --threads 4 -r $REGION $VCF1K |\
    bcftools view --threads 4 --samples-file $IDS --force-samples |\
    bcftools norm --threads 4 -m - - |\
    bcftools view --threads 4 --min-ac 2:minor - |\
    bcftools view --threads 4 -e INFO/AC==INFO/AN - |\
    bcftools annotate --threads 4 -x FORMAT -O z -o $VCFOUT -

plink --vcf $VCFOUT --r2 square spaces --real-ref-alleles --threads 4 --out ${PREFIX}
rm ${PREFIX}.nosex
