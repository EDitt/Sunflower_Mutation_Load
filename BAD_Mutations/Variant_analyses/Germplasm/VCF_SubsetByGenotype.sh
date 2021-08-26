#!/bin/bash

#SBATCH --job-name=VCF_subsetGenotype
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=12:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

module load BCFtools/1.10.2-GCC-8.3.0

# To run, need to define the following variable:
# 'group' - what group of samples to subset on
# 'alt_vcf' - the vcf subset containing sites where the alternate allele is deleterious
# 'ref_vcf' - the vcf subset containing sites where the reference allele is deleterious
# 'OUTPUT_DIR' - the directory to store results

mkdir -p ${OUTPUT_DIR}/${group}_sets

SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv

# 1. list of samples in that group
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)

# 2a. Subset using that list of samples to sites where at least 1 sample has the alternate allele:
bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="alt") > 0' -Oz -o ${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz # needs to be indexed for bcftools

NumAlt=$(bcftools view -Ov ${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz | grep -v "#" | wc -l)
echo "The group ${group} has at least one deleterious alternate allele at ${NumAlt} sites"

# 2b. Subset using that list of samples to sites where at least 1 sample has the reference allele:
#### note: the GT="ref" here only refers to homozygous sites, which is why we also have to add '| GT="het"
####        this is not the case for GT="alt" which refers to homozygous OR heterozygous implicitly
bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="ref"| GT="het") > 0' -Oz -o ${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz

NumRef=$(bcftools view -Ov ${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz | grep -v "#" | wc -l)
echo "The group ${group} has at least one deleterious reference allele at ${NumRef} sites"

# 3. Concatenate the VCFs
bcftools concat -a \
${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz \
${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz \
-Oz -o ${OUTPUT_DIR}/${group}_AllDel.vcf.gz

tabix -p vcf ${OUTPUT_DIR}/${group}_AllDel.vcf.gz

# 4. Clean up
# rm ${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz
# rm ${OUTPUT_DIR}/${group}_sets/${group}_AltDel.vcf.gz.tbi
# rm ${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz
# rm ${OUTPUT_DIR}/${group}_sets/${group}_RefDel.vcf.gz.tbi


