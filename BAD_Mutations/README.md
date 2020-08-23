# Predicting deleterious variants using BAD_Mutations: https://github.com/MorrellLAB/BAD_Mutations

## Variant Effect Predictor (VeP)

#### Prepare input files
Sort gff3 file and index with tabix
```bash
module load SAMtools/1.10-GCC-8.3.0
cd /scratch/eld72413/Ha412HOv2.0
grep -v "#" Ha412HOv2.0-20181130.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Ha412HOv2.0.gff3.gz
tabix -p gff Ha412HOv2.0.gff3.gz
```

### Install
```bash

cd /usr/local/singularity-images/

singularity exec ./ensembl-vep.simg vep [options]

singularity exec ./ensembl-vep.simg vep [options]
```

Make directory for local machine to store cache data for VEP
```bash
/scratch/eld72413/NSFproj/VEP_cacheFiles

docker run -t -i -v /scratch/eld72413/NSFproj/VEP_cacheFiles:/opt/vep/.vep ensemblorg/ensembl-vep
```
```bash
git clone https://github.com/Ensembl/ensembl-vep
cd ensembl-vep

module load Perl/5.26.0-GCCcore-6.4.0
export PATH=$PATH:/home/eld72413/perl5/bin/Archive/Zip
export PATH=$PATH:/home/eld72413/perl5/bin/Archive
git clone https://github.com/redhotpenguin/perl-Archive-Zip.git
perl Makefile.PL PREFIX=/home/eld72413/perl5/bin LIB=/home/eld72413/perl5/bin
make
make test
make install

# set up "cache files"
perl INSTALL.pl --CACHEDIR /scratch/eld72413/NSFproj/VEP_cacheFiles --NO_HTSLIB
#failed without --no_HTSLIB flag

# you do not have Archive:Zip installed

```
You do not have 'Archive::Zip' installed - Please install it as soon as possible
Non-zero exit status: 255
  Parse errors: No plan found in TAP output
Files=41, Tests=1680, 112 wallclock secs ( 0.18 usr  0.09 sys + 103.55 cusr  6.66 csys = 110.48 CPU)
Result: FAIL
Failed 1/41 test programs. 0/1680 subtests failed.


 - unpacking ./Bio/tmp/biodbhts.zip to ./Bio/tmp/
./Bio/tmp/Bio-DB-HTS-2.11 - moving files to ./biodbhts
 - making Bio::DB:HTS
Can't locate Module/Build.pm in @INC (@INC contains: ./Bio /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at Build.PL line 20.
BEGIN failed--compilation aborted at Build.PL line 20.
ERROR: Shared Bio::DB:HTS library not found

```bash

singularity exec ./ensembl-vep.simg vep [options]

singularity exec ./ensembl-vep.simg vep -i examples/homo_sapiens_GRCh38.vcf --cache

./vep -i examples/homo_sapiens_GRCh38.vcf --cache
```

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