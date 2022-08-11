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
#SBATCH --array=2-287


module load R/4.0.0-foss-2019b

module load BCFtools/1.13-GCC-8.3.0

module load BEDTools/2.30.0-GCC-8.3.0

# calculates number of SNPs in a specific genomic region (i.e. ROH) for a specific genotype, across variant classes

# variables that need to be defined
### SampleFile - a tab-delimited key with a list of samples with the first two columns corresponding to Plink's FID and IID, and the 3rd column corresponding to the gentype name in the VCF
### outdir - the output directory to save results
### ROH - the output of plink with each ROH for each individual
### vcf - the vcf file to use (can be a subset with only dSNP and sSNP)

mkdir -p ${outdir}/intermediates

mkdir -p ${outdir}/intermediates/GenotypeFiles

mkdir -p ${outdir}/intermediates/SNPs_ROH

#1.) get the sample names:
declare -a sample_array=($(awk '{print $3}' "${SampleFile}"))
declare -a plink_name_array=($(awk '{print $2}' "${SampleFile}"))

Sample="${sample_array[${SLURM_ARRAY_TASK_ID}]}"
Plink_Sample="${plink_name_array[${SLURM_ARRAY_TASK_ID}]}"

#2.) make a bedfile for the runs of homozygosity, including roh length
awk -v var="$Plink_Sample" '{OFS="\t"}; {if ($2==var) {print $4 "\t" $7-1 "\t" $8 "\t" $9}}' $ROH | \
sort -k 1,1 -k2,2n > ${outdir}/intermediates/ROH_${Sample}.bed

#3.) get variant alleles in those regions for each sample
bcftools view ${vcf} -Ou -s ${Sample} -R ${outdir}/intermediates/ROH_${Sample}.bed | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' |\
awk '{if ($5>0) {print $0}}' > ${outdir}/intermediates/ROH_${Sample}_SNPstats.txt

#4.) convert to derived and ancestral
Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Scripts/Variant_Table.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_Positions.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"${outdir}/intermediates/GenotypeFiles" \
"${Sample}" \
"${outdir}/intermediates/ROH_${Sample}_SNPstats.txt"

#5.) # put all SNPs that are homozygous, derived into a .bed file & combine with ROH lengths
### columns: chromosome, position, variant type, roh_start_position (0-based), roh_end_position, roh length (Kb)
awk '{if ($12==1) {print $1"\t"$2-1"\t"$2"\t"$13}}' ${outdir}/intermediates/GenotypeFiles/${Sample}_SNP_info.txt |\
sort -k 1,1 -k2,2n |\
bedtools intersect -a - -b ${outdir}/intermediates/ROH_${Sample}.bed -wao |\
awk '{OFS="\t"}; {print $1,$3,$4,$6,$7,$8}' > ${outdir}/intermediates/SNPs_ROH/${Sample}_SNP_ROH.txt

