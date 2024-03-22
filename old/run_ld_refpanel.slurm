#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=8:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=LD1KGP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/LD1KGP-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/LD1KGP-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

cd $SLURM_SUBMIT_DIR

CPUS=$SLURM_CPUS_PER_TASK

IDSALL=./data/kgp_eur_eas.txt
VCF1K=/reference_databases/1000G_VCF/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
PREFIX=./data/chr1_eur+eas
VCFOUT=${PREFIX}.vcf.gz

#bcftools view --threads $CPUS -r 1:200000000-210000000 $VCF1K |\
#    bcftools view --threads $CPUS --include '%ID=@./data/gwas_rsid.txt' - |\
#    bcftools view --threads $CPUS --samples-file $IDSALL --force-samples |\
#    bcftools view --threads $CPUS --genotype ^miss --min-ac 1:minor - |\
#    bcftools view --threads $CPUS -e INFO/AC==INFO/AN - |\
#    bcftools norm --threads $CPUS -m -both -O z -o $VCFOUT -
#
#plink --vcf $VCFOUT --r square spaces --real-ref-alleles --out ${PREFIX}
#rm ${PREFIX}.nosex

IDSEUR=./data/kgp_eur.txt
PREFIXEUR=./data/chr1_eur
VCFOUTEUR=${PREFIXEUR}.vcf.gz

bcftools view --threads $CPUS -r 1:200000000-210000000 $VCF1K |\
    bcftools view --threads $CPUS --include '%ID=@./data/gwas_rsid.txt' - |\
    bcftools view --threads $CPUS --samples-file $IDSEUR --force-samples |\
    bcftools view --threads $CPUS --genotype ^miss --min-ac 1:minor - |\
    bcftools view --threads $CPUS -e INFO/AC==INFO/AN - |\
    bcftools norm --threads $CPUS -m -both -O z -o $VCFOUTEUR -

plink --vcf $VCFOUTEUR --r square spaces --real-ref-alleles --out ${PREFIXEUR}
rm ${PREFIXEUR}.nosex

IDSEAS=./data/kgp_eas.txt
PREFIXEAS=./data/chr1_eas
VCFOUTEAS=${PREFIXEAS}.vcf.gz

bcftools view --threads $CPUS -r 1:200000000-210000000 $VCF1K |\
    bcftools view --threads $CPUS --include '%ID=@./data/gwas_rsid.txt' - |\
    bcftools view --threads $CPUS --samples-file $IDSEAS --force-samples |\
    bcftools view --threads $CPUS --genotype ^miss --min-ac 1:minor - |\
    bcftools view --threads $CPUS -e INFO/AC==INFO/AN - |\
    bcftools norm --threads $CPUS -m -both -O z -o $VCFOUTEAS -

plink --vcf $VCFOUTEAS --r square spaces --real-ref-alleles --out ${PREFIXEAS}
rm ${PREFIXEAS}.nosex

