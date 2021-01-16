#!/bin/bash

#SBATCH --job-name=Vcftools_freq
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL


VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/VeP/SAM_Sunflower_Subset.vcf.qz"
OUTPUTDIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/VeP"

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

vcftools --gzvcf $VCF \
--freq \
--out ${OUTPUTDIR}/SAM_SNPs_SUBSETFINAL
