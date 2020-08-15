#PBS -S /bin/bash
#PBS -q batch
#PBS -N SRA_Seq_Download
#PBS -l nodes=1:ppn=6
#PBS -l walltime=24:00:00
#PBS -l mem=20gb
#PBS -M dittmare@gmail.com
#PBS -m abe
#PBS -t 1-7

###This script downloads sequence files from the SRA
#As of 2.9.1 fasterq-dump can be used

export PATH=/home/eld72413/apps/sra/sratoolkit.2.9.6-centos_linux64/bin:${PATH}
#module load SRA-Toolkit/2.9.1-centos_linux64

#List of SRR#'s to download. Include full filepath
SRALIST="/scratch/eld72413/NSFproj/ancestralseqs/Debilis/Debilis_SRA_list.txt"

OUTPUTDIR="/scratch/eld72413/NSFproj/ancestralseqs/Debilis/RawSeqs"

TEMPDIR="/scratch/eld72413/NSFproj/ancestralseqs/Temp" #so temporary files are not saved in home directory (not enough space!)

SAMPLE=$(sed -n ${PBS_ARRAYID}p $SRALIST) #finds the sample #

fasterq-dump $SAMPLE \
-O $OUTPUTDIR \
-t $TEMPDIR \
-f \
--skip-technical
