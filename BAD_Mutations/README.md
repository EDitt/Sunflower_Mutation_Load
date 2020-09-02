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
Ran VeP.sh on cultivated and wild H. annuus VCF files using VeP.sh script

Filtered output
```bash
module load VEP/95.0-foss-2018b-Perl-5.28.0
INPUT=/scratch/eld72413/NSFproj/VEP/fullsam_remappedHa412HO_all.txt
filter_vep -i ${INPUT} -o Missense/fullsam_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT} -o Synon/fullsam_synon.txt -filter "Consequence is synonymous_variant"

INPUT2=/scratch/eld72413/NSFproj/VEP/wild_env_remappedHa412HO_all.txt
filter_vep -i ${INPUT2} -o Missense/wildenv_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT2} -o Synon/wildenv_synon.txt -filter "Consequence is synonymous_variant"
```

Filter VCF files
```bash
# list of positions
grep -v "#" fullsam_missense.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > fullsam_missense_positions.txt
grep -v "#" fullsam_synon.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > fullsam_synon_positions.txt
SAMVCF=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.fullsam.90.bi.remappedHa412HO_reheader.vcf.gz

module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1

vcftools --gzvcf $SAMVCF --freq --positions fullsam_missense_positions.txt --out SAM_MISSENSE
# After filtering, kept 119,332 out of a possible 2155376 Sites
vcftools --gzvcf $SAMVCF --freq --positions fullsam_synon_positions.txt --out SAM_SYNON
# After filtering, kept 184,668 out of a possible 2155376 Sites

grep -v "#" wildenv_missense.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > wildenv_missense_positions.txt
grep -v "#" wildenv_synon.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > wildenv_synon_positions.txt

WILDANN=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.env.90.bi.remappedHa412HO_reheader.vcf.gz

vcftools --gzvcf $WILDANN --freq --positions wildenv_missense_positions.txt --out WILD_MISSENSE
# After filtering, kept 306,904 out of a possible 4,882,321 Sites

vcftools --gzvcf $WILDANN --freq --positions wildenv_synon_positions.txt --out WILD_SYNON
# After filtering, kept 468,294 out of a possible 4,882,321 Sites

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