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

# LD pruning with Plink

module load PLINK/1.9b_5-x86_64

# variables that need to be defined:
### Output_Dir: the output directory for which to save the SNP lists
### File_Prefix: prefix for the .ped and .map files
### Window_Size: the size of window within which to prune SNPs (in Kb)
### Step_Size: the number of SNPs for which to shift the window at each step
### Rsquared: the R^2 threshold above which SNPs within those windows will be pruned


mkdir -p ${Output_Dir}/PrunedLists

plink --file ${File_Prefix} \
--indep-pairwise ${Window_Size} kb ${Step_Size} ${Rsquared} \
--allow-extra-chr \
--set-missing-var-ids @:# \
--out "${Output_Dir}"/PrunedLists/Plink_"${Rquared}"
