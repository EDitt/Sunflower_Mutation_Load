#!/bin/bash

#SBATCH --job-name=Filter_Variants_forGWAS
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=72:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -e
set -o pipefail

#### Filtering for GWAS

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

module load BCFtools/1.13-GCC-8.3.0

module load picard/2.16.0-Java-1.8.0_144
PICARD_JAR=/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar

# Variables to specify on command line:
# INPUT_VCF
# GEN_FASTA
# OUTPUT_DIR
# TEMP_DIR
# OUT_PREFIX

mkdir -p ${OUTPUT_DIR}/Intermediates

# Step 0: Check if there is a sequence dictionary. If not, make one
if [[ -f "${GEN_FASTA%.fasta}.dict" ]]; then
	echo "Sequence dictionary found, proceeding to step 1"
else
	echo "No sequence dictionary found. Creating one."
	java -jar ${PICARD_JAR} CreateSequenceDictionary \ 
      R="${GEN_FASTA}" \
      O="${GEN_FASTA%.fasta}.dict"
fi

# Step 1: Removing sites that failed variant recalibrator

if [[ -f "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_recalibrated_pass_sites.vcf.idx" ]]; then
	echo "Recalibrated pass sites already selected, proceeding to step 2"
else
	echo "Selecting Pass Sites from Variant Recalibrator"
	gatk SelectVariants \
		-R ${GEN_FASTA} \
		-V ${INPUT_VCF} \
		--exclude-filtered true \
		--create-output-variant-index true \
		--select-type-to-include SNP \
		--exclude-non-variants true \
		--tmp-dir ${TEMP_DIR} \
		-O ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_recalibrated_pass_sites.vcf

	num_sites1=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_recalibrated_pass_sites.vcf" | wc -l)

	echo "After filtering out sites that failed variant recalibrator, there are ${num_sites1} sites left"
fi

# Step 2: Removing sites with more than 0.1 missing genotypes
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf ]]; then
	echo "Sites with missing genotypes already filtered out, proceeding to step 3"
else
	echo "Removing sites with more than 0.1 missing genotypes"
	bcftools filter -e 'F_MISSING > 0.1' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_recalibrated_pass_sites.vcf > ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf
	num_sites2=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf" | wc -l)
	echo "After filtering out sites with too many low quality or missing variants, there are ${num_sites2} sites left"
fi

# Step 3: Filtering out sites with more than 0.1 heterozygous genotypes
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf ]]; then
	echo "Sites with more than 0.1 heterozygous genotypes already filtered out, proceeding to step 4"
else
	echo "Filtering out sites with more than 0.1 heterozygous genotypes"
	bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.1' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf -o ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf
	num_sites3=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf" | wc -l)
	echo "After filtering out sites with more than 0.1 heterozygous sites, there are ${num_sites3} sites remaining"
fi

# Step 4: Filtering out sites with < 0.01 minor allele frequency
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MAF_Filtered.vcf ]]; then
	echo "Sites with more <1% minor allele frequency already filtered out, proceeding to step 5"
else
	echo "Filtering out sites with <1% minor allele frequency"
	bcftools filter -i 'MAF > 0.01' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf -o ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MAF_Filtered.vcf
	num_sites4=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MAF_Filtered.vcf" | wc -l)
	echo "After filtering out sites with <1% minor allele frequency, there are ${num_sites4} sites remaining"
fi

# Step 5: Select only biallelic sites
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_BIALLELIC.vcf.gz ]]; then
	echo "Biallelic sites already selected"
else
	echo "Selecting biallelic sites"
	bcftools view -m2 -M2 -v snps --threads 4 ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_MAF_Filtered.vcf --output-type z --output-file ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_BIALLELIC.vcf.gz
fi

# Step 6: Index VCF
echo "Indexing VCF"
bcftools index ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_BIALLELIC.vcf.gz

# Step 7: Normalize VCF
echo "Normalizing VCF"
bcftools norm ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_BIALLELIC.vcf.gz \
--check-ref ws \
--fasta-ref ${GEN_FASTA} \
--threads 4 \
--output ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC_NORM.vcf.gz \
--output-type z

# Step 8: Indexing normalized VCF
echo "Indexing final VCF"
bcftools index ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC_NORM.vcf.gz

# Step 9: Get stats on final vcf:
bcftools stats --threads 4 ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC_NORM.vcf.gz > ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC_NORM_stats.txt
