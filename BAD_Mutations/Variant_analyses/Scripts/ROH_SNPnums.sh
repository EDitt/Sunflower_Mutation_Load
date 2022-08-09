#!/bin/bash

#SBATCH --job-name=ROH_SNP_info
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-287


# calculates number of SNPs in a specific genomic region (i.e. ROH) for a specific genotype, across variant classes

# variables that need to be defined
### SampleFile - a tab-delimited key with a list of samples with the first two columns corresponding to Plink's FID and IID, and the 3rd column corresponding to the gentype name in the VCF
### outdir - the output directory to save results
### ROH - the output of plink with each ROH for each individual
### vcf - the vcf file to use (can be a subset with only dSNP and sSNP)

mkdir -p ${outdir}/intermediates

mkdir -p ${outdir}/intermediates/GenotypeFiles

#1.) get the sample names:
declare -a sample_array=($(awk '{print $3}' "${SampleFile}"))
declare -a plink_name_array=($(awk '{print $2}' "${SampleFile}"))

Sample="${sample_array[${SLURM_ARRAY_TASK_ID}]}"
Plink_Sample="${plink_name_array[${SLURM_ARRAY_TASK_ID}]}"

#2.) make a bedfile for the runs of homozygosity
awk -v var="$Plink_Sample" '{OFS="\t"}; {if ($2==var) {print $4 "\t" $7-1 "\t" $8}}' $ROH | \
sort -k 1,1 -k2,2n > ${outdir}/intermediates/ROH_${Sample}.bed

#3.) get info for each allele in sample
bcftools view ${vcf} -Ou -s ${Sample} -R ${outdir}/intermediates/ROH_${Sample}.bed | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' |\
awk '{if ($5>0) {print $0}}' > ${outdir}/intermediates/ROH_${Sample}_SNPstats.txt

#4.) convert to derived and ancestral
Rscript "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/Variant_Table.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_Positions.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"${outdir}/intermediates/GenotypeFiles" \
"${Sample}" \
"${outdir}/intermediates/ROH_${Sample}_SNPstats.txt"

