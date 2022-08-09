#!/bin/bash

#SBATCH --job-name=Runs_homozyg
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# detecting runs of homozygosity with Plink

module load PLINK/1.9b_5-x86_64

plink --file ${File_Prefix} \
--homozyg \
--homozyg-group \
--allow-extra-chr \
--homozyg-window-kb 5000 \
--homozyg-window-snp 50 \
--homozyg-window-het 1 \
--homozyg-window-missing 5 \
--homozyg-window-threshold 0.05 \
--homozyg-snp 100 \
--homozyg-kb 1000 \
--homozyg-density 50 \
--homozyg-gap 1000 \
--out ${Output_Dir}/${Output_Prefix}
