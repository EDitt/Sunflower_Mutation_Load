#!/bin/bash

#SBATCH --job-name=Filter_VariantsExcessHet
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL


# this script will filter out variants with high (>5) values of Excess Heterozygosity

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

GEN_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter4_010621/Sunflower_SAM_SNP_Calling_HETFiltered.vcf
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter5_011121
TEMP_DIR="/scratch/eld72413/Tmp"

# First will need to index the file
gatk --java-options "-Xmx22g" IndexFeatureFile \
-F $INPUT_VCF

gatk --java-options "-Xmx22g" SelectVariants \
-R ${GEN_FASTA} \
-V ${INPUT_VCF} \
-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_EXCESSHETFiltered.vcf \
--selectExpressions "ExcessHet < 5.0" \
--tmp-dir ${TEMP_DIR}
