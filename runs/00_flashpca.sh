#$ -cwd
#$ -l h_rt=01:00:00
#$ -l h_vmem=4G
#$ -N GenomicPCA
#$ -pe smp 16
#$ -o runs/logs/genomic_pca/genomic_pca.out
#$ -e runs/logs/genomic_pca/genomic_pca.err

PCADIR=/home/home01/scrb/01_repos/GWAS_pipeline/data/transforms/GenomicPCA
BFILE=${PCADIR}/genotypes/only_cmr/ukb_cal_all_chrs_v2_CMR_GBR_indiv
FLASHPCA=/home/home01/scrb/nobackup/software/flashpca_x86-64

$FLASHPCA \
  --bfile $BFILE 
  --ndim 20 -v \
  --outpc ${PCADIR}/genomic_pcs.txt

	 
