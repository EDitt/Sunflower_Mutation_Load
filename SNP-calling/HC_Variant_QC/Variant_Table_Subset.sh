#PBS -S /bin/bash
#PBS -q batch
#PBS -N VCF_Stats
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -l mem=20gb

#PBS -M dittmare@gmail.com
#PBS -m abe

set -o pipefail


RAW_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Sunflower_SAM_SNP_Calling_raw_variants.vcf
HC_SUBSET=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Sunflower_SAM_SNP_Calling_filtered.vcf
OUTPUTDIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/QC_filter

module load BCFtools/1.9-foss-2016b
module load vcflib/20180410

#First subset raw VCF file
echo "Randomly subsetting raw VCF file to 1%"
bcftools view "${RAW_VCF}" | vcfrandomsample -p 65 -r 0.01 > "${OUTPUTDIR}/Raw_variants_subsample.vcf"

# intersection of raw variants with hc_subset (get the overlap)
# need to compress vcf files for the bcftools command...
#bcftools isec -n=2 "${OUTPUTDIR}/Raw_variants_subsample.vcf" "${HC_SUBSET}" > "${OUTPUTDIR}/"
module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1
echo "Comparing sites between high-confidence subset and raw VCF sub-sample"
vcftools --vcf "${OUTPUTDIR}/Raw_variants_subsample.vcf" --diff "${HC_SUBSET}" --diff-sites --out "${OUTPUTDIR}/HC_Raw_sites"

# turn both into tables (include whether SNPs or Indels)
# get complement in R (positions not shared)