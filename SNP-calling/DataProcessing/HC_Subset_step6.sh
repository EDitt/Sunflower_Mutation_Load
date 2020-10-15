#PBS -S /bin/bash
#PBS -q batch
#PBS -N HC_Subset_step6
#PBS -l nodes=1:ppn=6
#PBS -l walltime=72:00:00
#PBS -l mem=50gb

#PBS -M dittmare@gmail.com
#PBS -m abe

module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1

vcfoutput="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Sunflower_SAM_SNP_Calling_filtered.vcf"
OUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset"

echo "Removing sites that aren't polymorphic."

vcftools --vcf "${vcfoutput}" --non-ref-ac 1 --recode --recode-INFO-all --out "${OUT_DIR}/Sunflower_SAM_SNP_Calling_high_confidence_subset"
mv "${OUT_DIR}/Sunflower_SAM_SNP_Calling_high_confidence_subset.recode.vcf" "${OUT_DIR}/Sunflower_SAM_SNP_Calling_high_confidence_subset.vcf" # Rename the output file

echo "Finished removing sites that aren't polymorphic."