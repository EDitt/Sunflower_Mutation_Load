#!/bin/bash

#SBATCH --job-name=GWAS_filters
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -e
set -o pipefail


module load BCFtools/1.10.2-GCC-8.3.0

mkdir -p ${OUTPUT_DIR}/Intermediates

# filter for missing data, minor allele frequency, and too many heterogygous individuals
### VARIABLES TO SPECIFY:
##### INPUT_VCF - VCF to filter
##### OUTPUT_DIR - Output directory
##### OUT_PREFIX - Prefix name for VCF files
##### 

# Step 1: Heterozygosity (<10%) & Minor Allele frequency (>1%)

bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.1 && MAF > 0.01' \
${INPUT_VCF} \
-o ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MissHet_Filtered.vcf

num_sites1=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MissHet_Filtered.vcf" | wc -l)

echo "After filtering out sites with MAF <1% and heterozygosity >10%, there are ${num_sites1} sites left"

# Step 2: Missingness (<10%)

bcftools filter -e 'F_MISSING > 0.1' \
${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MissHet_Filtered.vcf > ${OUTPUT_DIR}/${OUT_PREFIX}_missing_filtered.vcf

num_sites2=$(grep -v "#" "${OUTPUT_DIR}/${OUT_PREFIX}_missing_filtered.vcf" | wc -l)

echo "After filtering out sites with >10% missing data, there are ${num_sites2} sites left"
