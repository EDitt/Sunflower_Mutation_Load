#!/bin/bash

#SBATCH --job-name=GzipVCF
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# a script to gzip and index a vcf file. Must specify file variable in command line:
#sbatch --export=file=<file> gzip_vcf.sh

module load SAMtools/1.10-iccifort-2019.5.281

bgzip -c --threads 4 $file > ${file}.gz
tabix -p vcf ${file}.gz
