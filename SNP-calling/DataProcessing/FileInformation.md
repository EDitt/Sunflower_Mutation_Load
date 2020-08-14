### Raw Sequence Data Information

Re-sequencing data from 288 samples lines

##### Raw sequence data Location:
Most sequences were available in /project/jmblab/sunflower_1_sequence_data

Sequences were separated into groups for processing. Copying by group was done using the following code:
```bash
OUTPUTDIR="/scratch/eld72413/SAM_seq/results2/Group3/RawSeqs"
INPUTDIR="/project/jmblab/sunflower_1_sequence_data"
LIST=/scratch/eld72413/SAM_seq/results2/Mapping_redo3.txt

GROUP=3 #or whichever group was being copied

for line in `awk -v var="$GROUP" '{if ($2 == var) {print $1}}' $LIST`; do
	cp -R $INPUTDIR/${line} $OUTPUTDIR
done
```

Some sequences (N=9) needed to be obtained from the SRA:  
Sequence data for the following samples were not available in the Sunflower_1_sequence_data folder:  
PPN267
PPN271
PPN274
PPN279
PPN281
PPN282
PPN283
PPN287
PPN288  
Sequence data for these samples was downloaded from the SRA using the SRA toolkit on 4/25/19.

Data were first retrieved using the `prefetch` command, and then converted to fastq files using the `fastq-dump` command with the "split-files" flag to obtain forward and reverse sequences.
Then, names were changed from SRA numbers to PPN numbers. Fastq files were then gzipped for temporary storage.

Each file contained only one forward and reverse sequence, except for PPN282 which contained 2 forward and 2 reverse sequences. For this sample, the forward and reverse sequences were merged using the cat command.

Some sequences (N=8) were in the folder: /project/jmblab/sunflower_1_sequence_data/South_Africa_seqs:
PPN012                  
PPN022                  
PPN035                  
PPN057                  
PPN060                  
PPN068                  
PPN167                  
PPN175                  
PPN227                  
PPN233                  
These samples had many associated sequence files (83-378)

# Check integrity of copied files:

##### Pre-processing
Many samples had more than 1 forward/reverse sequence file that needed to be concatenated before beginning sequence handling
