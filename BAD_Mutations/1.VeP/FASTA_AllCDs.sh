#!/bin/bash

#SBATCH --job-name=FASTA_cds
#SBATCH --partition=batch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o FastaCDs.out
#SBATCH -e FastaCDs.err

# Makes 1 FASTA file with all CDs

module load gffread/0.11.6-GCCcore-8.3.0

GFF3=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3
FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files

gffread $GFF3 -g $FASTA -x ${OUTPUTDIR}/All_CDs.fasta