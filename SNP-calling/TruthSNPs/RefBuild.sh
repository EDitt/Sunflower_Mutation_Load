#PBS -S /bin/bash
#PBS -q batch
#PBS -N bowtie_BuildRef
#PBS -l mem=48gb,nodes=1:ppn=4,walltime=06:00:00
#PBS -m abe
#PBS -M dittmare@gmail.com

set -e
set -o pipefail

module load Bowtie2/2.3.5.1-GCC-8.3.0

# User provided input arguments
REF_DIR=/scratch/eld72413/Ha412HOv2.0
REF_FILENAME=Ha412HOv2.0-20181130.fasta

# Go into reference dir
cd ${REF_DIR}

# Generate database prefix
PREFIX=$(basename ${REF_FILENAME} .fasta)

# Create bowtie2 index database
bowtie2-build ${REF_DIR}/${REF_FILENAME} ${PREFIX}

# Check the content of the database
bowtie2-inspect -n ${PREFIX}