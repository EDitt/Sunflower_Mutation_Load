#!/bin/bash

#SBATCH --job-name=Pi_calc
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-16

module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

# variables that need to be defined:
### GenomeFile: a text file containing a list of chromosomes being analyzed (anything other than the first column will be ignored)
### VCF: the vcf file to calculate Fst
### GenoList: a text file list of genotypes to include
### OUT_PREFIX: The output directory and prefix to be used for output file(s)

declare -a chrom_array=($(awk '{print $1}' "${GenomeFile}"))

CHROM="${chrom_array[${SLURM_ARRAY_TASK_ID}]}"


vcftools --gzvcf ${VCF} \
--keep ${GenoList} \
--chr ${CHROM} \
--window-pi 1000000 \
--window-pi-step 1000000 \
--out ${OUT_PREFIX}_${CHROM}

vcftools --gzvcf ${VCF} \
--keep ${GenoList} \
--chr ${CHROM} \
--TajimaD 1000000 \
--out ${OUT_PREFIX}_${CHROM}
