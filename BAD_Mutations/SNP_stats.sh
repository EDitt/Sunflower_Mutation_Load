#!/bin/bash

#SBATCH --job-name=SNP_stats
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=12:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

module load BCFtools/1.10.2-GCC-8.3.0

SNP_set="/scratch/bp26123/UBC_snps/Annuus.tranche90.snp.remappedHa412.vcf.gz"
OUTPUT_DIR="/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset_Raw"
#tabix $SNP_set

bcftools stats --threads 4 $SNP_set > ${OUTPUT_DIR}/UBC_Dataset_Stats.txt

