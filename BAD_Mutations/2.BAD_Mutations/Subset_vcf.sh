#!/bin/bash

#SBATCH --job-name=Subset_vcf
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

module load BCFtools/1.10.2-GCC-8.3.0
bcftools view -R ${positions} ${vcf} > ${outputdir}/${name}.vcf
