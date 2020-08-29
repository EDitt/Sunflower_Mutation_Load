#PBS -S /bin/bash
#PBS -q batch
#PBS -N Concatenate_SRA
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -l mem=2gb
#PBS -t 1-132

#PBS -M dittmare@gmail.com
#PBS -m abe

LIST=/home/eld72413/DelMut/Sunflower_Mutation_Load/Outgroups/WildAnnuusOnePerPop
INPUTDIR=/scratch/eld72413/NSFproj/ancestralseqs/Annuus/Raw
OUTPUTDIR=/scratch/eld72413/NSFproj/ancestralseqs/Annuus/Raw/Merged

Sample=$(sed -n ${PBS_ARRAYID}p | awk '{print $2}')

NumFiles=$(find -maxdepth 1 -name "${Sample}*1.fastq" | wc -l)

if [[ "$NumFiles" -eq 0 ]]
then
echo "No $Sample Files in $INPUTDIR"
else
echo "Concatenating $NumFiles Forward $Sample"
cat `ls ${INPUTDIR}/${Sample}*1.fastq` > ${OUTPUTDIR}/"$Sample"_R1.fastq
echo "Concatenating $NumFiles Reverse $Sample"
cat `ls ${INPUTDIR}/${Sample}*2.fastq` > ${OUTPUTDIR}/"$Sample"_R2.fastq
fi
