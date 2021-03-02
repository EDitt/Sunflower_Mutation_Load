#!/bin/bash

#SBATCH --job-name=Filter_Variants_AllSteps
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=48:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -e
set -o pipefail

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

module load BCFtools/1.10.2-GCC-8.3.0

GEN_FASTA="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All"
#OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test"
TEMP_DIR="/scratch/eld72413/Tmp"
OUT_PREFIX="Sunflower_SAM_SNP"

### Cut-offs
GQ=6
MinDP=3
MaxDP=50
QUAL=40.0
het_prop=0.2 # I get errors when I try to use the variable in the bcftools filter statement
MaxExcessHet=5.0

mkdir -p ${OUTPUT_DIR}/Intermediates

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

# Step 2: Removing low quality genotypes (min/max DP and min GQ) and sites with more than 0.2 missing genotypes

if [[ -f "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoField.vcf.idx" ]]; then
	echo "Filter labels already added to low quality genotypes"
else
	echo "Adding Filter labels to low quality genotypes. Cut-offs are as follows: MinGQ is ${GQ}; MinDP is ${MinDP}; MaxDP is ${MaxDP}"
	gatk --java-options "-Xmx22g" VariantFiltration \
		-R ${GEN_FASTA} \
		-V ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_recalibrated_pass_sites.vcf \
		-O ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoField.vcf \
		--genotype-filter-name "lowGQ" --genotype-filter-expression "GQ < ${GQ}"  \
		--genotype-filter-name "highDP" --genotype-filter-expression "DP > ${MaxDP}" \
		--genotype-filter-name "lowDP" --genotype-filter-expression "DP < ${MinDP}" \
		--tmp-dir ${TEMP_DIR}

	echo "Done adding filter labels to low quality genotypes"
fi

if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoFieldFiltered.vcf.idx ]]; then
	echo "Filtered genotypes already set to no call, proceeding to step 3"
else
	echo "Removing sites with too much missing data, then setting filtered genotypes to no call"
	gatk --java-options "-Xmx22g" SelectVariants \
	-R ${GEN_FASTA} \
	-V ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoField.vcf \
	-O ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoFieldFiltered.vcf \
	--set-filtered-gt-to-nocall true \
	--max-fraction-filtered-genotypes 0.2 \
	--max-nocall-fraction 0.2 \
	--exclude-non-variants true \
	--remove-unused-alternates true \
	--tmp-dir ${TEMP_DIR}

	num_sites2=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoFieldFiltered.vcf" | wc -l)
	echo "After filtering out sites with too much missing variants, there are ${num_sites2} sites left"
fi

### NOTE: It was discovered that the order of filtering by gatk in the previous step was that it removed nocall-fraction genotypes *before* it set filtered genotypes to no-call,
# necessitating an intermediate step here to filter out sites with too many filtered genotypes - used bcftools (see bcftools_missingfilter.sh)
# Step 3 below was added later when this was discovered:

# Step 3: Filter out missing genotypes (GATK sets filtered genotypes to nocall *after* filtering for missing data)
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf ]]; then
	echo "Sites with too many low quality genotypes (set to missing by GATK in previous step) already filtered out, proceeding to step 4"
else
	echo "Removing sites with more than 0.2 low quality genotypes"
	bcftools filter -e 'F_MISSING > 0.2' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_GenoFieldFiltered.vcf > ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf
	num_sites3=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf" | wc -l)
	echo "After filtering out sites with too many low quality or missing variants, there are ${num_sites3} sites left"
fi

# Step 4: Filtering out sites with too many heterozygous genotypes
if [[ -f ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf ]]; then
	echo "Sites with more than ${het_prop} heterozygous genotypes already filtered out, proceeding to step 5"
else
	echo "Filtering out sites with more than ${het_prop} heterozygous genotypes"
	bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.2' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf -o ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf
	num_sites4=$(grep -v "#" "${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf" | wc -l)
	echo "After filtering out sites with more than ${het_prop} heterozygous sites, there are ${num_sites4} sites remaining"
fi

# Step 5: Filtering out sites with low QUAL values or high ExcessHet values
if [[ -f ${OUTPUT_DIR}/${OUT_PREFIX}_FINALFilter.vcf.idx ]]; then
	echo "Sites with QUAL values below ${QUAL} and/or ExcessHet values above ${MaxExcessHet} already filtered out, proceeding to step 6"
else
	echo "Removing sites with QUAL values below ${QUAL} and/or ExcessHet values above ${MaxExcessHet}"
# first index file
	gatk --java-options "-Xmx2g" IndexFeatureFile \
		-F ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf

	gatk --java-options "-Xmx22g" SelectVariants \
		-R ${GEN_FASTA} \
		-V ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf \
		-O ${OUTPUT_DIR}/${OUT_PREFIX}_FINALFilter.vcf \
		--selectExpressions "ExcessHet < ${MaxExcessHet}" \
		--selectExpressions "QUAL > ${QUAL}" \
		--exclude-non-variants true \
		--tmp-dir ${TEMP_DIR}

	num_sites5=$(grep -v "#" "${OUTPUT_DIR}/${OUT_PREFIX}_FINALFilter.vcf" | wc -l)
	echo "After filtering for QUAL and ExcessHet annotations, there are ${num_sites5} sites remaining"
fi

# Step 6: Select only biallelic sites
if [[ -f ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC.vcf ]]; then
	echo "Biallelic sites already selected"
else
	echo "Selecting biallelic sites"
	bcftools view -m2 -M2 -v snps --threads 4 ${OUTPUT_DIR}/${OUT_PREFIX}_FINALFilter.vcf --output-type v --output-file ${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC.vcf
	num_sites6=$(grep -v "#" "${OUTPUT_DIR}/${OUT_PREFIX}_BIALLELIC.vcf" | wc -l)
	echo "After selecting only biallelic sites, there are ${num_sites6} variants"
fi
