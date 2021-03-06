#!/bin/bash

#SBATCH --job-name=VeP
#SBATCH --partition=highmem_p

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL


#Compressed_VCF=/scratch/eld72413/NSFproj/PublishedSNPs/Edited/fullsam.90.remappedHa412HO_norm_biallelic.vcf.gz
VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz"
OUTPUTDIR="/scratch/eld72413/SAM_seq/VeP"
OUTPUTPREFIX=SAM_SNP_Final_BiallelicNorm

#    Variant sets should be either 'deletions', 'insertions', or 'snps'
VARIANT_SET=all

module load VEP/101.0-foss-2019b-Perl-5.30.0
# cd /usr/local/singularity-images/
# singularity exec ./ensembl-vep.simg vep \

vep \
    -i ${VCF} \
    --gff /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0.gff3.gz \
    --fasta  /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta \
    --species  helianthus_annuus \
    --check_svs \
    --verbose \
    --format vcf \
    --warning_file ${OUTPUTDIR}/${OUTPUTPREFIX}_WARN.txt \
    -o ${OUTPUTDIR}/${OUTPUTPREFIX} \
    --buffer_size 1000000 \
    --fork 4 \
    --force_overwrite

## can use flag --buffer_size [number] to decreae run time. (# of variants that are read into memory). Default=5000
## flag --no_intergenic does not include intergenic consequences in the output