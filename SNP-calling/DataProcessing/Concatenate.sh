#PBS -S /bin/bash
#PBS -q batch
#PBS -N Concatenate
#PBS -l nodes=1:ppn=6
#PBS -l walltime=80:00:00
#PBS -l mem=20gb
#PBS -t 1-50

#PBS -M dittmare@gmail.com
#PBS -m abe

INPUTDIR=/scratch/eld72413/SAM_seq/results2/Group4/RawSeqs
OUTPUTDIR=/scratch/eld72413/SAM_seq/results2/Group4/MergedSeqs

DIR=$(find $INPUTDIR -name "PPN*" | sed -n ${PBS_ARRAYID}p)

if [[ -d "$DIR" ]]; then
	name=$(basename ${DIR})
	for f1 in `ls "$DIR"/*R1.fastq.gz`; do
		if [[ -f "$f1" ]]; then
			f2=${f1%%1.fastq.gz}"2.fastq.gz"
			echo "Concatenating $name R1"
			zcat $f1 >> ${OUTPUTDIR}/"$name"_R1.fastq.gz
			echo "Concatenating $name R2"
			zcat $f2 >> ${OUTPUTDIR}/"$name"_R2.fastq.gz
		fi
	done
fi