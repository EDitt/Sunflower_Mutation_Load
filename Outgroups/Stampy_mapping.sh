#PBS -S /bin/bash
#PBS -q batch
#PBS -N Map_Outgroup
#PBS -l nodes=1:ppn=24
#PBS -l walltime=24:00:00
#PBS -l mem=100gb
#PBS -M dittmare@gmail.com
#PBS -m abe
#PBS -t 1-7


### This script will map reference sequence files to a genome using Stampy

set -o pipefail


module load Stampy/1.0.31-foss-2016b-Python-2.7.14

###### User provided arguments

#	Directory containing Trimmed sequences to map
SEQ_DIR='/scratch/eld72413/NSFproj/ancestralseqs/Debilis/Adapter_Trimming'

#   Where do our output files go?
OUT_DIR='/scratch/eld72413/NSFproj/ancestralseqs/Debilis/Mapped/05Perc'

#	Naming scheme for forward and reverse reads
FORWARDNAME="_Forward_ScytheTrimmed.fastq.gz"

REVERSENAME="_Reverse_ScytheTrimmed.fastq.gz"

#   What is our per site substitution rate?
DIVERGENCE='0.05'

#   Reference prefix must match .stidx files
REF_PREFIX='Ha412HOv2'

#   What directory contains our reference files?
REF_DIR='/scratch/eld72413/NSFproj/ancestralseqs/GenomeFiles'

#### End of User provided arguments

#	Finds forward read
f1=$(find ${SEQ_DIR} -maxdepth 1 -name "*${FORWARDNAME}" | sed -n ${PBS_ARRAYID}p)

#	Finds Reverse read
f2=${f1%%$FORWARDNAME}"$REVERSENAME"

#   Takes our forward read and generates sample name for output file
SAMPLE_NAME=$(basename ${f1%%$FORWARDNAME}"")

#   Stampy
#       -g is the genome index file, PREFIX.stidx
#       -h is the genome hash file, PREFIX.sthash
#		-t is the number of threaads to use
#       --substitutionrate is the expected fraction of Poisson-distributed substitutions (default is 0.001)
#       -f is the output file format (SAM is default)
#       -M reads to map (needs to be on the last line or else forward/reverse reads need to be separated by a comma)
#       -o is our output file

cd ${REF_DIR}
stampy.py -g "${REF_PREFIX}" \
       -h "${REF_PREFIX}" \
       -t 8 \
       --substitutionrate="${DIVERGENCE}" \
       -f sam \
       -o "${OUT_DIR}/${SAMPLE_NAME}_${DIVERGENCE}.sam" \
       -M "${f1}","${f2}"