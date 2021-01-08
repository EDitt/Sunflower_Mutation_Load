#!/bin/bash

#SBATCH --job-name=Filter_VariantsHets
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

INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter3_123120/Sunflower_SAM_SNP_Calling_QUALFiltered.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter4_010621"

bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.2' $INPUT_VCF -o ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_HETFiltered.vcf
