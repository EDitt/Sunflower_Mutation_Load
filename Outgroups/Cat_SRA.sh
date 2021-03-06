#PBS -S /bin/bash
#PBS -q batch
#PBS -N Concatenate_SRA
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -l mem=2gb
#PBS -t 1-16

#PBS -M dittmare@gmail.com
#PBS -m abe

LIST=/scratch/eld72413/NSFproj/ancestralseqs/Landrace/RawSeqs/names2
INPUTDIR=/scratch/eld72413/NSFproj/ancestralseqs/Landrace/RawSeqs
OUTPUTDIR=/scratch/eld72413/NSFproj/ancestralseqs/Landrace/RawSeqs/Merged

Sample=$(sed -n ${PBS_ARRAYID}p $LIST)

NumFiles=$(find $INPUTDIR -maxdepth 1 -name "${Sample}_*1.fastq" | wc -l)

F_Files=$(ls ${INPUTDIR}/${Sample}_*1.fastq)
R_Files=$(ls ${INPUTDIR}/${Sample}_*2.fastq)

if [[ "$NumFiles" -eq 0 ]]; then
	echo "No $Sample Files in $INPUTDIR"
else
	echo "Concatenating $NumFiles Forward $Sample"
	echo "Forward Files: $F_Files"
	cat $F_Files > ${OUTPUTDIR}/"$Sample"_R1.fastq
	echo "Concatenating $NumFiles Reverse $Sample"
	echo "Reverse Files: $R_Files"
	cat $R_Files > ${OUTPUTDIR}/"$Sample"_R2.fastq
fi
