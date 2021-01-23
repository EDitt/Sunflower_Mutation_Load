#!/bin/bash

#SBATCH --job-name=makeFASTA_cds
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o MakeFastas.sh.%A_%a.out
#SBATCH -e MakeFastas.sh.%A_%a.err

#SBATCH --array=1-57


module load gffread/0.11.6-GCCcore-8.3.0
# gffread will make a FASTA file out of the CDS that falls within a specified genomic region

GFF3=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3
FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
OUTPUTDIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Biallelic/VeP/FASTAs

List_to_Process=$(sed -n ${SLURM_ARRAY_TASK_ID}p $List_of_lists)
echo "Processing List ${List_to_Process}"

while read line; do
	coord=$(grep "${line}" $GFF3 | awk '{if ($3=="mRNA") {print $1":"$4".."$5}}')
	num=$(for i in ${coord}; do echo $i; done | wc -l)
	if [[ $num > 1 ]]; then
		echo "ERROR: There are ${num} lines matching for ${line}"
	else
		echo "Processing FASTA for ${line}"
		gffread $GFF3 -g $FASTA -r ${coord} -x ${OUTPUTDIR}/${line#mRNA:}.fasta
		sed -i 's/>mRNA:/>/' ${OUTPUTDIR}/${line#mRNA:}.fasta
	fi
done < $List_to_Process
