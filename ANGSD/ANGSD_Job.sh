#!/bin/bash

#SBATCH --job-name=ANGSD
#SBATCH --partition=highmem_p

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --tmp=10g

set -o pipefail

module load GSL/2.6-iccifort-2019.5.281
#module load SAMtools/1.10-GCC-8.3.0
#module load gnuplot/5.2.2-foss-2018a
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC

cd /home/eld72413/ANGSD_DEV/new3_afterMaintenance/angsd-wrapper

./angsd-wrapper ${WRAPPER} ${CONFIG}
