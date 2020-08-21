# Predicting deleterious variants using BAD_Mutations: https://github.com/MorrellLAB/BAD_Mutations

## Variant Effect Predictor (VeP)

### Install
```bash
git clone https://github.com/Ensembl/ensembl-vep
cd ensembl-vep

# set up "cache files"
perl INSTALL.pl

```
 - unpacking ./Bio/tmp/biodbhts.zip to ./Bio/tmp/
./Bio/tmp/Bio-DB-HTS-2.11 - moving files to ./biodbhts
 - making Bio::DB:HTS
Can't locate Module/Build.pm in @INC (@INC contains: ./Bio /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at Build.PL line 20.
BEGIN failed--compilation aborted at Build.PL line 20.
ERROR: Shared Bio::DB:HTS library not found

## Prep

I cloned this repository, and then switched to `dev` branch:
```bash
git checkout -t origin/dev
```

Modules needed:
```bash
module load Biopython/1.74-foss-2018a-Python-2.7.14
#module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+/2.10.0
module load HyPhy/2.5.15-gompi-2019b
module load PASTA/1.8.2-foss-2016b-Python-2.7.14

# PASTA, HyPhy, Clustal-omega, fasttree
```

## Setup

Variables
```bash
OUTPUTDIR="/scratch/eld72413/NSFproj/BADMutations"
```

```bash
# to show available species databases
python /home/eld72413/DelMut/BAD_Mutations/BAD_Mutations.py setup --list-species
# helianthus not in this list

python /home/eld72413/DelMut/BAD_Mutations/BAD_Mutations.py -v DEBUG setup -c $OUTPUTDIR -b $OUTPUTDIR -t 'Hannuus' -d /home/eld72413/apps
# getting errors
```
I tried creating my own config file based off of Sample_Config.txt

```bash
TEST_CONFIG=/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Test_Config081720
```

## Fetch

```bash
python /home/eld72413/DelMut/BAD_Mutations/BAD_Mutations.py -v DEBUG fetch -c $TEST_CONFIG 
```