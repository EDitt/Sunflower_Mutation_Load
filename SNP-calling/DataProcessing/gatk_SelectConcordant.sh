#!/bin/bash

#SBATCH --job-name=Concordance
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

GEN_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset_Raw/SAM_lines/UBC_Dataset_SAMlines.recode.vcf"
TRUTH_VCF="/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf"
OUTPUT_DIR="/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset_Raw/SAM_lines/Consensus"
TEMP_DIR="/scratch/eld72413/Tmp"

# this will output variants that are also in the truth set. Then I can see what their annotation values are to see if my filtering was too stringent

gatk --java-options "-Xmx22g" SelectVariants \
-R ${GEN_FASTA} \
-V ${INPUT_VCF} \
-O ${OUTPUT_DIR}/UBC_SAM_TruthConcordant.vcf \
--concordance ${TRUTH_VCF} \
--tmp-dir ${TEMP_DIR}
