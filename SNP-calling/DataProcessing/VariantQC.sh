#!/bin/bash

#SBATCH --job-name=Variant_Table
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

INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed


echo "Using GATK's variants to table function"
gatk VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F TYPE -F QUAL -F FILTER -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F TRANSITION -F HET \
     -O "${OUTPUT_DIR}/Variants.table" \
     --show-filtered
echo "done"

#    -GF GQ -GF DP \
