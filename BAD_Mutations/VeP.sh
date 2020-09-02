#PBS -S /bin/bash
#PBS -q batch
#PBS -N VEP
#PBS -l nodes=1:ppn=1
#PBS -l walltime=18:00:00
#PBS -l mem=22gb
#PBS -m abe
#PBS -M dittmare@gmail.com

#module load HTSlib/1.9-foss-2018b
#module load perl/modules.centos7.5.26.1
#module load BCFtools/1.6-foss-2016b

#     The VCF
#Compressed_VCF=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.fullsam.90.bi.remappedHa412HO_reheader.vcf.gz
#Compressed_VCF=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.env.90.bi.remappedHa412HO_reheader.vcf.gz
Compressed_VCF=/scratch/eld72413/NSFproj/VEP/Synon/SAM_SYNON.recode.vcf.gz
OUTPUTDIR=/scratch/eld72413/NSFproj/VEP/Synon
OUTPUTPREFIX=sam_synonymous

#    Variant sets should be either 'deletions', 'insertions', or 'snps'
VARIANT_SET=all

module load VEP/95.0-foss-2018b-Perl-5.28.0
# cd /usr/local/singularity-images/
# singularity exec ./ensembl-vep.simg vep \

vep \
    -i ${Compressed_VCF} \
    --gff /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0.gff3.gz \
    --fasta  /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta \
    --species  helianthus_annuus \
    --total_length \
    --check_svs \
    --verbose \
    --format vcf \
    --warning_file ${OUTPUTDIR}/${OUTPUTPREFIX}_WARN.txt \
    -o ${OUTPUTDIR}/${OUTPUTPREFIX}

