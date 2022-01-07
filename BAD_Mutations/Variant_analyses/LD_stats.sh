#!/bin/bash

#SBATCH --job-name=LD_stats
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

# calculate R^2 and D' statistics for pairs of SNPs

module load PLINK/1.9b_5-x86_64

# variables that need to be defined:
### GenomeFile: a text file containing a list of chromosomes being analyzed (anything other than the first column will be ignored)
### File_Prefix: prefix for the .ped and .map files
### Kb_Window_Size: the size of window for pairs of SNPs for which to calculate r^2 (default 1000)
### NumVariant_Windows: the maximum number of variants between two SNPs allowed for r^2 calculations to be performed on them (default 10)
### MinR2_Window: the minimum R2 number, below which will not be reported between variants (default 0.2)

declare -a chrom_array=($(awk '{print $1}' "${GenomeFile}"))

CHROM="${chrom_array[${SLURM_ARRAY_TASK_ID}]}"

plink --r2 inter-chr dprime \
--file ${File_Prefix} \
--allow-extra-chr \
--chr ${CHROM} \
--ld-window-kb ${Kb_Window_Size} \
--ld-window ${NumVariant_Windows} \
--ld-window-r2 ${MinR2_Window} \
--out ${CHROM}.r0
