#!/bin/bash -l
#PBS -q batch
#PBS -N VEP
#PBS -l nodes=1:ppn=1
#PBS -l walltime=18:00:00
#PBS -l mem=4gb
#PBS -m abe
#PBS -M dittmare@gmail.com

​
module load HTSlib/1.9-foss-2018b
#module load perl/modules.centos7.5.26.1
module load BCFtools/1.6-foss-2016b

#     The original VCF
VCF=.vcf
Compressed_VCF=.vcf.gz
​
#    Compress the VCF file
#bgzip --index ${VCF}
#bcftools view ${VCF} --output-file ${Compressed_VCF} -output-type z
​
#while [ ! -f ${Compressed_VCF} ]; do sleep 5; done
​
#bcftools index --force ${Compressed_VCF}
​
#    Variant sets should be either 'deletions', 'insertions', or 'snps'
VARIANT_SET=all
​
module avail VEP/95.0-foss-2018b-Perl-5.28.0
cd 
​
/home/morrellp/shared/Software/ensembl-vep-release-97.3/vep \
    -i ${Compressed_VCF} \
    --gff /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3 \
    --fasta  /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta \
    --species  \
    --total_length \
    --check_svs \
    --verbose \
    --format vcf \
    --warning_file Morex_Mutants_${VARIANT_SET}.txt \
    -o Morex_Mutants_${VARIANT_SET}.txt