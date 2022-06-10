#!/bin/bash

#SBATCH --job-name=GenotypeStats
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-7


# outputs files with SNP counts for all variant classes

### Variables: 
# Table: (SNP_INFO), i.e. starting table (can be all snps or subset that are polarized)
# outputdir: output directory
# vcf: full VCF

#######################

module load BCFtools/1.13-GCC-8.3.0

mkdir -p ${outputdir}/intermediates

declare -a annotation_array=($(awk '{print $13}' $Table | sort -u))

VariantClass="${annotation_array[${SLURM_ARRAY_TASK_ID}]}"


# subset alt and ref derived

awk -v var="$VariantClass" 'BEGIN{FS=OFS="\t"}; {if ($3 == $11 && $13 == var) {print $1,$2}}' ${Table} > ${outputdir}/intermediates/${VariantClass}_Alt_derivedPositions.txt

awk -v var="$VariantClass" 'BEGIN{FS=OFS="\t"}; {if ($4 == $11 && $13 == var) {print $1,$2}}' ${Table} > ${outputdir}/intermediates/${VariantClass}_Ref_derivedPositions.txt

# output sample counts for all samples

### note: the bcftools stats "-s -" flag means to include all samples
bcftools view -Oz ${vcf} -R ${outputdir}/intermediates/${VariantClass}_Alt_derivedPositions.txt | \
bcftools stats -s - | \
grep "PSC" > ${outputdir}/intermediates/${VariantClass}_Alt_derivedCounts.txt

bcftools view -Oz ${vcf} -R ${outputdir}/intermediates/${VariantClass}_Ref_derivedPositions.txt | \
bcftools stats -s - | \
grep "PSC" > ${outputdir}/intermediates/${VariantClass}_Ref_derivedCounts.txt




