#PBS -S /bin/bash
#PBS -q batch
#PBS -N bowtie_Align
#PBS -l mem=48gb,nodes=1:ppn=6,walltime=06:00:00
#PBS -m abe
#PBS -M dittmare@gmail.com

set -e
set -o pipefail

module load Bowtie2/2.3.5.1-GCC-8.3.0

# This script does a quick alignment as a check for the reference genome
# Following what Li did for morex v1: https://github.com/lilei1/9k_BOPA_SNP/blob/master/script/commandlines

# User provided input arguments
FASTA_FILE=/scratch/eld72413/SNParray/ContextualSeqs_UniqMap.fasta
DB_NAME=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130
OUT_PREFIX=Contextual_UniqMap_HA412v2_bowtie2
OUT_DIR=/scratch/eld72413/SNParray/AllUniqSNPs

# Check if our dir exists, if not make it
#mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with bowtie2
bowtie2 -f ${FASTA_FILE} -x ${DB_NAME} -S ${OUT_PREFIX}.sam