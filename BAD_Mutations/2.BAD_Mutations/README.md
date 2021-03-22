# Using BAD_Mutations on UMN cluster

## Testing 1 region

### Set up software environment
```bash
ssh mesabi #need to do this to load modules or conda environments
module load python3/3.6.3_anaconda5.0.1
module load parallel #will need this later when parallelizing
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations
```

Made a directory for Sunflower data:
`/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/`
working in `Tests` sub-directory to test code on 2 regions

### Create config file
The database of CDS files has already been created in: `/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes`

```bash
cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations
# how is sunflower named in database?
./BAD_Mutations.py setup --list-species

./BAD_Mutations.py setup \
	-b /panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes \
	-t 'Hannuus' \
	-e 0.05 \
	-c /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/config.txt
```
Result:
```bash
===2021-01-23 00:07:04,755 - Setup_Env===
INFO	Wrote configuration into /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/config.txt
```

### Align

```bash
srun -N 1 -n 1 -c 1 --mem=8gb -t 12:00:00 -p interactive --pty bash #639487
tmux new -s align
module load python3/3.6.3_anaconda5.0.1
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations
./BAD_Mutations.py -v DEBUG align \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/config.txt \
    -f /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161.fasta \
    -o /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161 2> /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Alignment.log
```
This took 10 min for 1 region

### Predict

```bash
tmux new -s predict
srun -N 1 -n 1 -c 1 --mem=8gb -t 12:00:00 -p interactive --pty bash #640374
module load python3/3.6.3_anaconda5.0.1
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations
./BAD_Mutations.py -v DEBUG predict \
	-c /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/config.txt \
	-f /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161.fasta \
	-a /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Ha412HOChr12g0573161_MSA.fasta \
	-r /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Ha412HOChr12g0573161.tree \
	-s /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161.subs \
	-o /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Prediction \
	1> /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Prediction/Predict.log
```

```bash
===2021-01-23 01:49:53,854 - LRT_Prediction===
DEBUG   Aligned Pos: 3486, 3489, 3492, 4659, 5103, 5307, 5322, 5325, 5532, 6309, 6318, 6831, 7014
===2021-01-23 01:49:53,855 - LRT_Prediction===
DEBUG   HyPhy input file: 
/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Ha412HOChr12g0573161_MSA.fasta
/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Ha412HOChr12g0573161.tree
/tmp/BAD_Mutations_HYPHY_Subs_pv8b6g8o.txt
Ha412HOChr12g0573161
===2021-01-23 01:49:53,856 - LRT_Prediction===
DEBUG   bash /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations/Shell_Scripts/Prediction.sh /home/morrellp/liux1299/.conda/envs/bad_mutations/bin/HYPHYMP /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations/Shell_Scripts/LRT.hyphy /tmp/BAD_Mutations_HYHPY_In_9bq3ai4t.txt /tmp/BAD_Mutations_HYPHY_Out_2u5sepgn.txt
===2021-01-23 02:14:37,134 - LRT_Prediction===
DEBUG   stdout:

===2021-01-23 02:14:37,135 - LRT_Prediction===
DEBUG   stderr:

Check messages.log details of this run.

===2021-01-23 02:14:37,155 - LRT_Predict===
INFO    Prediction in /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Prediction/Ha412HOChr12g0573161_Predictions.txt
```

### Compile
```bash
cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations
./BAD_Mutations.py compile \
	--pred-dir /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161/Prediction \
	--long-subs /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Tests/Ha412HOChr12g0573161.subs
```

Error message: (wiki should be edited to reflect the required commands)
```
usage: BAD_Mutations.py compile [-h] --pred-dir PRED_DIR --long-subs LONG_SUBS
BAD_Mutations.py compile: error: the following arguments are required: --pred-dir/-P, --long-subs/-S
```

## Entire dataset

Uploaded .tar.gz archive to MSI (`~shared/Projects/Sunflower`) as well as the Bad Mutaions summary file
(See README.md from 1.VeP for how files were prepared)

Extracted
```bash
cd /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower
tar -xf Bad_mutations.tar.gz
```

## Split CDs into separate FASTA files
Using Chaochih's script (split_cds_fasta.py from Barley_Mutated.git repository) to split up the CDs into individual files

```bash
module load python3_ML/3.7.1_anaconda
./split_cds_fasta.py /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/All_CDs.fasta /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/BadMutationsFASTAS
```
Loading fasta file...
Splitting fasta file...
Done.

Rename so extension is .fasta and not .fa & take out the "mRNA:" prefix
Also check number of sequences and number of characters
```bash
sbatch --export=FASTA_DIR='/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/BadMutationsFASTAS' FASTA_check.sh # 1592462

grep "ERROR:" FastaCheck.out | wc -l # 530 (not divisible by 3)
```

### Create config file
The database of CDS files has already been created in: `/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes`

Conda software environment has been setup: `/home/morrellp/liux1299/.conda/envs/bad_mutations`

```bash
ssh mesabi # no longer need to do this because no login anymore on MSI
module load python3/3.6.3_anaconda5.0.1
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations
./BAD_Mutations.py setup \
    -b /panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes \
    -t 'Hannuus' \
    -e 0.05 \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/config.txt
```

Result:
```bash
===2021-01-27 13:52:22,612 - Setup_Env===
INFO	Wrote configuration into /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/config.txt
```

### Align

Chaochih Liu created scripts for parallelizing the alignment and prediction steps across many regions. See: https://github.com/MorrellLAB/Barley_Mutated/tree/master/02_analysis/bad_mutations

Chaochih created a list of the individual FASTA files for each CD and then split that into lists containing 500 sequence records for each file, then created a list of the lists

Create a list of all fasta files, break into sublists of 500 each, make list of these lists
```bash
FASTA_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/BadMutationsFASTAs"
find ${FASTA_DIR} -name "Ha*.fasta" | sort -V > /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/fasta_lists/All_Fasta.txt #50838

split -l 500 --numeric-suffixes All_Fasta.txt Hannuus_cds_list- --suffix-length=3 --additional-suffix=.txt #102 files

find $(pwd -P) -name "*list-*.txt" | sort -V > all_cds_Hannuus_list_of_lists.txt
```

Used Chaochih's script: `bad_mut_align.sh`. I only changed the User-defined input variables
```bash
sbatch --array=0-101 bad_mut_align.sh
```
sbatch: Setting account: morrellp
Submitted batch job 1604694


check number of files
```bash
cd /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/MSA_output
find -name "*.tree" | wc -l # 50492 
find -name "*_MSA.fasta" | wc -l # 50492 (was 41665 previously with FASTA issues)

Input=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/MSA_output

find $Input -maxdepth 1 -type d -name "*list-*" | while read dir; do
cd $dir
filenum=$(find -name "*.tree" | wc -l)
if [[ $filenum -ne 500 ]]; then
echo "$dir has $filenum files"
fi
done
# none of them have 500 files (most missing 1-9)
```

How many are scaffold sequence?
```bash
# in the FASTA directory
find -maxdepth 1 -name "Ha412HOChr00c*" | wc -l # 34

# in the MSA output directory
find -name "Ha412HOChr00c*tree" | wc -l #32


```
9173 regions are missing. It turns out this is because the way I was using gffread to make FASTA files would sometimes put multiple sequences in the same file
One sequence that failed: Ha412HOChr00c00023g0859521.fasta
I need to edit some of the FASTA files (see MakeFastas.sh and README in 1.VeP folder) because some FASTA files had > 1 sequence
I re-ran these again after re-making the FASTA files

Now only 346 are missing
-manually BLAST FASTA files


##
(on previous SNP set- JobID: 740887): array 71 had errors ("more than 1 record found in handle"). Need to check others for this error. (still running after 9 hours) Ha412HOChr12g0550871.fasta
- jobs 740887_107 and 740887_96 stayed on pending ("launch failed requeued held"). Deleted these jobs


### Predict

```bash
# make list of MSA output directories
cd /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/MSA_output
find $(pwd -P) -mindepth 1 -maxdepth 1 -type d -name "Hannuus*" | sort -V > MSA_output_dir_list.txt

# edit name of substitution file
cp ./Sunflower/SAM_SNP_BadMut_Summary ./SAM_SNP_BadMut_Summary_edit
SUB=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/SAM_SNP_BadMut_Summary_edit
sed -i 's/mRNA://' $SUB

# make file of transcripts (all are primary)
awk '{print $1}' ${SUB} | sort -u > /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Transcript_names.txt

```