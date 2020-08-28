#PBS -S /bin/bash
#PBS -q batch
#PBS -N ANGSD
#PBS -l nodes=1:ppn=6
#PBS -l walltime=24:00:00
#PBS -l mem=22gb

#PBS -M dittmare@gmail.com
#PBS -m abe

set -o pipefail

module load GSL/2.6-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load gnuplot/5.2.2-foss-2018a

WRAPPER=""
CONFIG=""

./angsd-wrapper ${WRAPPER} ${CONFIG}
