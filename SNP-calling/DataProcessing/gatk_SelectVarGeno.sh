#!/bin/bash

#SBATCH --job-name=Select_Variants
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

GEN_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter2_122828/Sunflower_SAM_SNP_Calling_GenoField.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter2_122828"
TEMP_DIR="/scratch/eld72413/Tmp"

# filtering: 
# --max-fraction-filtered-genotypes: maximum fraction of samples filtered at the genotype level
# --max-nocall-fraction: maximum fraction of samples with no-call genotypes
# non-variant sites (--exclude-non-variants )
# I want to set filtered genotypes to no-call- I *think* the --exclude-filtered flag only applies to INFO field,
# but --set-filtered-gt-to-nocall applies to genotype field

gatk --java-options "-Xmx22g" SelectVariants \
-R ${GEN_FASTA} \
-V ${INPUT_VCF} \
-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_GenoFieldFiltered.vcf \
--set-filtered-gt-to-nocall true \
--max-fraction-filtered-genotypes 0.2 \
--max-nocall-fraction 0.2 \
--exclude-non-variants true \
--remove-unused-alternates true \
--tmp-dir ${TEMP_DIR}
