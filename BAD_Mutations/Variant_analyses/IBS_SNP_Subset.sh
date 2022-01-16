#!/bin/bash

#SBATCH --job-name=IBS_forSNP_subset
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

# calculates pairwise IBS from a list of SNPs

module load PLINK/1.9b_5-x86_64

mkdir -p $outputdir/Intermediates

# variables that need to be defined:
## 1.) Input_prefix - the prefix for plink-formatted vcf (including directory)
## 2.) SNP_List - List of SNPs to use ('Chromosome:Position'). One per line. (Include directory)
## 3.) Window_Size
## 4.) Step_Size
## 5.) Rsquared
## 6.) outputdir - output directory for files
## 7.) output_prefix - output prefix for final files


#### 1.) LD prune SNP subset

plink --file  ${Input_prefix} \
--extract ${SNP_List} \
--indep-pairwise ${Window_Size} kb ${Step_Size} ${Rsquared} \
--allow-extra-chr \
--out ${outputdir}/Intermediates/TEMP_Plink_${Rsquared}

Num_removed=$(wc -l ${outputdir}/Intermediates/TEMP_Plink_${Rsquared}.prune.out)
Num_kept=$(wc -l ${outputdir}/Intermediates/TEMP_Plink_${Rsquared}.prune.in)

echo "Pruning SNPs with R^2 more than ${Rsquared} in ${Window_Size} kb windows will remove ${Num_removed} variants, keeping ${Num_kept} variants"

### combine pruned list with input list
comm -12 --check-order <(sort ${outputdir}/Intermediates/TEMP_Plink_${Rsquared}.prune.in) <(sort ${SNP_List}) > ${outputdir}/Intermediates/TEMP_SNPsToKeep.txt

Num_SNPs=$(wc -l ${outputdir}/Intermediates/TEMP_SNPsToKeep.txt)

echo "Pairwise IBS will be calculated using ${Num_SNPs} SNPs"

#### 2.) Calculate pairwise IBS with the subset of SNPs
 
plink --genome \
--file ${Input_prefix} \
--extract ${outputdir}/Intermediates/TEMP_SNPsToKeep.txt \
--allow-extra-chr \
--out ${outputdir}/Intermediates/TEMP_${output_prefix}

#### 3.) Clean up file for use (fix names that were altered in Plink conversion)

awk '{print $2, $4, $12}' ${outputdir}/Intermediates/TEMP_${output_prefix}.genome > ${outputdir}/${output_prefix}_IBS.txt

sed -i 's/PPN285/Hopi_PPN285/g' ${outputdir}/${output_prefix}_IBS.txt
sed -i 's/33/SF_33/g' ${outputdir}/${output_prefix}_IBS.txt
sed -i 's/PPN136/NMS373_PPN136/g' ${outputdir}/${output_prefix}_IBS.txt
sed -i 's/PPN251/RHA415-4_PPN251/g' ${outputdir}/${output_prefix}_IBS.txt
sed -i 's/531071/PI_531071/g' ${outputdir}/${output_prefix}_IBS.txt



