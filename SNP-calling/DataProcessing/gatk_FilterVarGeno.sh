#!/bin/bash

#SBATCH --job-name=Filter_Variants
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=12:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

REF_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120/Sunflower_SAM_SNP_Calling_snps.filtered.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter2_122828"

gatk VariantFiltration \
	-R ${REF_FASTA} \
	-V ${INPUT_VCF} \
	-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_GenoField.vcf \
	--genotype-filter-name "below6GQ" --genotype-filter-expression "GQ < 6"  \
	--genotype-filter-name "above50DP" --genotype-filter-expression "DP > 50" \
	--tmp-dir ${TEMP_DIR}
