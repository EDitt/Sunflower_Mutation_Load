#!/bin/bash

#SBATCH --job-name=VCFtoPlink
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

module load PLINK/1.9b_5-x86_64

plink --vcf $INPUT_VCF \
--recode \
--allow-extra-chr \
--real-ref-alleles \
--keep-allele-order \
--out $OUT_PREFIX
