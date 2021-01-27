# Using BAD_Mutations on UMN cluster

Started with test regions

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

### Entire dataset

Uploaded .tar.gz archive to MSI (`~shared/Projects/Sunflower`)

Extracted
```bash
cd /panfs/roc/groups/9/morrellp/shared/Projects/Sunflower
tar -xf Bad_mutations.tar.gz

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
