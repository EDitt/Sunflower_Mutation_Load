#!/bin/bash

#SBATCH --job-name=Filter_VariantsMinDP
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=24:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# filtering out genotypes with a dp below 3

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

REF_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921"
TEMP_DIR="/scratch/eld72413/Tmp"

gatk --java-options "-Xmx22g" VariantFiltration \
	-R ${REF_FASTA} \
	-V ${INPUT_VCF} \
	-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_GenoDPField.vcf \
	--genotype-filter-name "below3DP" --genotype-filter-expression "DP < 3" \
	--tmp-dir ${TEMP_DIR}

gatk --java-options "-Xmx22g" SelectVariants \
-R ${REF_FASTA} \
-V ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_GenoDPField.vcf \
-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_DP_min3Filtered.vcf \
--set-filtered-gt-to-nocall true \
--max-fraction-filtered-genotypes 0.2 \
--max-nocall-fraction 0.2 \
--exclude-non-variants true \
--remove-unused-alternates true \
--tmp-dir ${TEMP_DIR}
