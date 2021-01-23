### Set up config file
```bash
module load Biopython/1.75-foss-2019b-Python-2.7.16
# how is sunflower named in database
BAD_Mut_Path="/home/eld72413/DelMut/BAD_Mutations"
python ${BAD_Mut_Path}/BAD_Mutations.py setup --list-species
# Hannuus

#tblastx, pasta, hyphy, clustalo, fasttree
module load BLAST+/2.10.1-gompi-2019b
module load HyPhy/2.5.15-gompi-2019b # this wasn't added to config
module load FastTree/2.1.11-GCCcore-8.3.0 # this wasn't added to config
# FastTree & HyPhy paths: 
/apps/mf/eb/all
# could add this manually to config?
module load pasta/1.8.2_conda

# no CLUSTALO module

python ${BAD_Mut_Path}/BAD_Mutations.py setup \
    -b /scratch/eld72413/NSFproj/CDs_database \
    -t 'Hannuus' \
    -e 0.05 \
    -c /scratch/eld72413/NSFproj/BADMutations/config.txt


Install hyphy and fasttree with conda?
```bash
module load Anaconda3/2020.02
conda install -c bioconda hyphy
conda config --add channels bioconda
```
^ did not work because defaulted to directory I couldn't download to

### Install PASTA
```bash
cd /home/eld72413/apps
git clone https://github.com/smirarab/pasta.git
git clone https://github.com/smirarab/sate-tools-linux.git
cd pasta/
python setup.py develop  --user
```

### Set up config file
```bash
module load Biopython/1.75-foss-2019b-Python-2.7.16
# how is sunflower named in database
BAD_Mut_Path="/home/eld72413/DelMut/BAD_Mutations"
python ${BAD_Mut_Path}/BAD_Mutations.py setup --list-species
# Hannuus

#tblastx, pasta, hyphy, clustalo, fasttree
module load BLAST+/2.10.1-gompi-2019b
module load HyPhy/2.5.15-gompi-2019b # this wasn't added to config
module load FastTree/2.1.11-GCCcore-8.3.0 # this wasn't added to config
# FastTree & HyPhy paths: 
/apps/mf/eb/all
# could add this manually to config?
module load pasta/1.8.2_conda

# no CLUSTALO module

python ${BAD_Mut_Path}/BAD_Mutations.py setup \
    -b /scratch/eld72413/NSFproj/CDs_database \
    -t 'Hannuus' \
    -e 0.05 \
    -c /scratch/eld72413/NSFproj/BADMutations/config.txt

```

### Fetch
```bash
tmux new -s fetch
srun --pty  -p inter_p  --mem=4G --nodes=1 --ntasks-per-node=4 --time=12:00:00 --job-name=qlogin /bin/bash -lm #1210297

python ${BAD_Mut_Path}/BAD_Mutations.py -v DEBUG fetch -c /scratch/eld72413/NSFproj/BADMutations/config.txt
```
