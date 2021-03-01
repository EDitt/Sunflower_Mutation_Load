#!/bin/bash

#SBATCH --job-name=MissingFilter
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=24:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL


module load BCFtools/1.10.2-GCC-8.3.0


#bcftools view -e 'F_MISSING > 0.3' ${vcf} > ${out_dir}/Variant_Filtering/test_new_filter/${out_prefix}_missing_filtered_rawVcf_v2.vcf works the same according to Chaochih
#bcftools filter -e 'F_PASS(GT="mis") > 0.3' ${vcf} > ${out_dir}/Variant_Filtering/test_new_filter/${out_prefix}_missing_filtered_rawVcf.vcf

bcftools filter -e 'F_MISSING > 0.2' ${vcf} > ${out_dir}/${out_prefix}_missing_filtered.vcf