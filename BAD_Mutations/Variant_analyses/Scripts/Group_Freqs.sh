#!/bin/bash

#SBATCH --job-name=GroupFreqs
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-5



## a script to get frequency distribution for different groups of germplasm

### Variables:
# SAM_INFO: spreadsheet with information about each line
# VCF: vcf to subset from
# outputdir: output directory

mkdir -p ${outputdir}/intermediates

declare -a group_array=($(awk 'BEGIN {FS=","}NR>1 {print $14}' $SAM_INFO | sort -u))

Group="${group_array[${SLURM_ARRAY_TASK_ID}]}"

genotypes=$(awk -v var="$Group" -F',' '{if ($14==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO | paste -sd,) 
bcftools view -Ou --samples ${genotypes} ${VCF} |\
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${outputdir}/intermediates/${Group}_SNP_freq.txt

Rscript "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/Variant_Table.R" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/FinalPositionFiles" \
"${outputdir}/intermediates/${Group}_SNP_freq.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"${outputdir}/${Group}_SNP_info.txt"

Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${outputdir}/${Group}_SNP_info.txt" \
"1.0" \
"0.05" \
"Derived_Freq" \
"${outputdir}/${Group}_DerivedFreq_Bins.txt"

