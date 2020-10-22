#PBS -S /bin/bash
#PBS -q batch
#PBS -N Filter_SNPs_GATK
#PBS -l nodes=1:ppn=6
#PBS -l walltime=36:00:00
#PBS -l mem=22gb

#PBS -M dittmare@gmail.com
#PBS -m abe

module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
GATK_JAR=/usr/local/apps/eb/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk

GEN_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120"
TEMP_DIR="/scratch/eld72413/Tmp"

# filtering: 
# indels (--select-type-to-include SNP), 
# sites that have been tagged as FAIL by ApplyVQSR (--exclude-filtered)
# non-variant sites (--exclude-non-variants )
gatk --java-options "-Xmx22g" SelectVariants \
-R ${GEN_FASTA} \
-V ${INPUT_VCF} \
-O ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_snps.filtered.vcf \
--select-type-to-include SNP \
--exclude-filtered true \
--exclude-non-variants true \
--tmp-dir ${TEMP_DIR}
