
#filename is mRNA:Ha412HOChr17g0858291.subs

# loop through list of files, take out .subs and mRNA/ prefix

GFF3_file

name=$(awk '{if ($3=="CDS") {print $0}}' Ha412HOv2.0-20181130.gff3 | grep mRNA:Ha412HOChr17g0858291 | awk '{print $1":"$4"-"$5}')

REF_FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
OUTPUTDIR=/scratch/eld72413/NSFproj/VEP/BADMutationsFastas

module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

gatk FastaReferenceMaker \
	-R "$REF_FASTA" \
	-O "$OUTPUTDIR/Ha412HOChr17g0858291.fasta" \
	-L "${name}"

# adds a "1 " to the sequence name- need to change
grep ">1 " Ha412HOChr17g0858291.fasta