#!/bin/bash -l

#SBATCH --job-name=ANGSD
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=22G
#SBATCH --time=6:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --chdir=/scratch/eld72413/Tmp

set -o pipefail

module load GSL/2.6-iccifort-2019.5.281
module load SAMtools/1.10-iccifort-2019.5.281
module load gnuplot/5.2.8-GCCcore-8.3.0

#module load SAMtools/1.10-GCC-8.3.0
#module load gnuplot/5.2.2-foss-2018a
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC

SLURM_TMPDIR="/scratch/eld72413/Tmp"

cd /home/eld72413/ANGSD_DEV/new4_recompileREDO/angsd-wrapper

./angsd-wrapper ${WRAPPER} ${CONFIG}

#./angsd-wrapper Inbreeding /scratch/eld72413/SAM_seq/ANGSD/Configuration_Files/Inbreeding_Coefficients_Config
