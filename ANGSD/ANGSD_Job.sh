#PBS -S /bin/bash
#PBS -q batch
#PBS -N ANGSD
#PBS -l nodes=1:ppn=6
#PBS -l walltime=48:00:00
#PBS -l mem=50gb

#PBS -M dittmare@gmail.com
#PBS -m abe

set -o pipefail

module load GSL/2.6-GCC-8.3.0
#module load SAMtools/1.10-GCC-8.3.0
#module load gnuplot/5.2.2-foss-2018a
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC

WRAPPER="SFS"
CONFIG="/scratch/eld72413/NSFproj/ANGSD_FILES/Chrom10/SFS_Chrom10_Config"

cd /home/eld72413/ANGSD_DEV/new2/angsd-wrapper

./angsd-wrapper ${WRAPPER} ${CONFIG}
