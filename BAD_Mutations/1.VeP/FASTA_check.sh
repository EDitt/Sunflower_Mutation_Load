#!/bin/bash

#SBATCH --job-name=CheckNRenameFASTAs
#SBATCH -p small,ram256g,ram1t,max

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=22G
#SBATCH --time=36:00:00
#SBATCH --export=None 
#SBATCH --mail-user=dittmare@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o FastaCheck.out
#SBATCH -e FastaCheck.err

# checks that FASTA files all only have 1 sequence per file, # of characters divisible by 3 and renames + takes out mRNA prefix

cd $FASTA_DIR

for file in `ls -1`; do
	num_char=$(grep -v "^>" $file | grep -Eo '[[:alnum:]]' | wc -l)
	if (( $num_char % 3 != 0 )); then
		echo "ERROR: $file is not Divisible by 3"
	fi
	num_seq=$(grep "^>" $file | wc -l)
	if [[ $num_seq -ne 1 ]]; then
		echo "ERROR: $num sequences in $file"
	else
		sed -i 's/>mRNA:/>/' ${file}
		echo "Re-naming $file to ${file#mRNA:}sta"
		mv $file "${file#mRNA:}sta"
	fi
done
