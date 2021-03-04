#!/bin/bash

#SBATCH --job-name=NormalizeVCF
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL



OUTPUTDIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All
INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_BIALLELIC.vcf.gz
FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta

module load BCFtools/1.10.2-GCC-8.3.0

bcftools norm ${INPUT_VCF} \
--check-ref s \
--fasta-ref ${FASTA} \
--threads 4 \
--output ${OUTPUTDIR}/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz \
--output-type z
