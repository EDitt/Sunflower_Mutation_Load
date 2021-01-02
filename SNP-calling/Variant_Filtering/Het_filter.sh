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

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

REF_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter3_123120/Sunflower_SAM_SNP_Calling_QUALFiltered.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121"
TEMP_DIR="/scratch/eld72413/Tmp"


echo "Marking heterozygote genotypes and removing previous filters..."
gatk --java-options "-Xmx22g" VariantFiltration \
	-R ${REF_FASTA} \
	-V ${INPUT_VCF} \
	-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_HetField.vcf \
	--genotype-filter-name "isHet" --genotype-filter-expression "isHet == 1"  \
	--invalidate-previous-filters true \
	--tmp-dir ${TEMP_DIR}
echo "Finished marking heterozygote genotypes and removing previous filters"

echo "Selecting sites with no more than 20% heterozygotes"
# filtering: 
gatk --java-options "-Xmx22g" SelectVariants \
-R ${REF_FASTA} \
-V ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_HetField.vcf \
-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_HetFieldFiltered.vcf \
--max-fraction-filtered-genotypes 0.20 \
--tmp-dir ${TEMP_DIR}
echo "Finished selecting sites with no more than 20% heterozygotes"
