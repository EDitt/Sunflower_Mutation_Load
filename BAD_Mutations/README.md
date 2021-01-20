# Predicting deleterious variants using BAD_Mutations: https://github.com/MorrellLAB/BAD_Mutations

# My SAM SNPs

## Prepare input files

This takes a really long time with full set of SNPs. To QC data, will first subset using my regions file with 1M random intervals of 200 bp regions
```bash
# for bcftools need a 3-column regions file
Bed="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed"
OUTPUTDIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/VeP"

awk '{print $1,$2,$3}' $Bed > ${OUTPUTDIR}/RegionsFile1M.bed

srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l
regions="${OUTPUTDIR}/RegionsFile1M.bed"
VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf.gz
module load BCFtools/1.10.2-GCC-8.3.0

bcftools filter -R ${Bed} ${VCF} -o ${OUTPUTDIR}/SAM_Sunflower_Subset.vcf.qz
bcftools stats ${OUTPUTDIR}/SAM_Sunflower_Subset.vcf.qz #3,264,946 variants
# number multi-allelic
awk '{if ($3 > 2) {print $0}}' SAM_SNPs_SUBSETFINAL.frq | wc -l #189,373

awk '{if ($3 == 2) {print $0}}' SAM_SNPs_SUBSETFINAL.frq > SAM_SNPs_SUBSETFINAL_Biallelic.frq
```

Normalize VCF (make sure allele matches reference)
Used `NormalizeVCF.sh` script

New index + stats on normalized vcf:
```bash
tmux new -s bcftools_biallelic
srun --pty  -p inter_p  --mem=4G --nodes=1 --ntasks-per-node=4 --time=6:00:00 --job-name=qlogin /bin/bash -l

module load BCFtools/1.10.2-GCC-8.3.0

bcftools index --threads 4 Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz

bcftools stats --threads 4 Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz > Sunflower_SAM_SNP_Calling_norm_biallelic_stats.txt

```



# UBC SNPs

## Prepare input files

```bash
#qsub -I -q s_interq -l walltime=12:00:00 -l nodes=1:ppn=4 -l mem=8gb

out_dir=/scratch/eld72413/NSFproj/PublishedSNPs/Edited
in_dir=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset
FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta

#module load SAMtools/1.9-foss-2016b-htslib
module load BCFtools/1.9-foss-2016b

bcftools norm ${in_dir}/Annuus.tranche90.snp.fullsam.90.bi.remappedHa412HO_reheader.vcf.gz \
--check-ref s \
--fasta-ref $FASTA \
--threads 4 \
--output ${out_dir}/fullsam.90.remappedHa412HO_norm.vcf.gz \
--output-type z
#Lines   total/split/realigned/skipped:  2155376/0/179/0
#REF/ALT total/modified/added:   2155376/340799/925019

bcftools index ${out_dir}/fullsam.90.remappedHa412HO_norm.vcf.gz

bcftools stats --threads 4 ${out_dir}/fullsam.90.remappedHa412HO_norm.vcf.gz > ${out_dir}/fullsam.90.remappedHa412HO_norm_norm.vcf.stats.txt 
# 2,155,376 SNPs

# Is this file already filtered for only biallelic SNPs?:
bcftools view -m 3 ${out_dir}/fullsam.90.remappedHa412HO_norm.vcf.gz | wc -l #950166 - no
# Filter for bi-allelic SNPs
bcftools view -m2 -M2 -v snps --threads 4 ${out_dir}/fullsam.90.remappedHa412HO_norm.vcf.gz --output-type z --output-file ${out_dir}/fullsam.90.remappedHa412HO_norm_biallelic.vcf.gz

# new index + stats
bcftools index ${out_dir}/fullsam.90.remappedHa412HO_norm_biallelic.vcf.gz
bcftools stats --threads 4 ${out_dir}/fullsam.90.remappedHa412HO_norm_biallelic.vcf.gz > ${out_dir}/fullsam.90.remappedHa412HO_norm_biallelic_stats.txt
# 1,230,357

```

#### Sort gff3 file and index with tabix
```bash
module load SAMtools/1.10-GCC-8.3.0
cd /scratch/eld72413/Ha412HOv2.0
grep -v "#" Ha412HOv2.0-20181130.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Ha412HOv2.0.gff3.gz
tabix -p gff Ha412HOv2.0.gff3.gz
```

## Variant Effect Predictor (VeP)

Ran VeP.sh on cultivated and wild H. annuus VCF files using VeP.sh script
From Peter's notes:
-    Using the `total_length` flag in VeP puts some extra information in the substitions file that BAD_Mutations does not like
-    Removing that manually to test!
-    This worked, so I've altered `ensembl_vep_Fagioli.sh` to run without that option. Need to regenerate substitutions and run the process again!

#### VeP filter (not needed for BAD_Mutations)

Filter output to look at distribution of different variant classes
```bash
module load VEP/95.0-foss-2018b-Perl-5.28.0
INPUT=/scratch/eld72413/NSFproj/VEP/fullsam_remappedHa412HO_all.txt
filter_vep -i ${INPUT} -o Missense/fullsam_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT} -o Synon/fullsam_synon.txt -filter "Consequence is synonymous_variant"

INPUT2=/scratch/eld72413/NSFproj/VEP/wild_env_remappedHa412HO_all.txt
filter_vep -i ${INPUT2} -o Missense/wildenv_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT2} -o Synon/wildenv_synon.txt -filter "Consequence is synonymous_variant"
```

#### Format VeP output for BAD_Mutations
```bash
VEP_OUTPUT=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/fullsam_remappedHa412HO_norm_biallelic

#gzip input file
gzip -c ${VEP_OUTPUT} > fullsam_remappedHa412HO_norm_biallelic.txt.gz

VEP_OUTPUT_GZIP=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/fullsam_remappedHa412HO_norm_biallelic.txt.gz
# output file needs to include all directory information
OUTPUTFILE=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/fullsam_remapped_norm_biallelic_BMsummary
# this will be a directory that contains substitution files for every transcript
OUTPUTDIR=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/sub_files

module load Biopython/1.74-foss-2018a-Python-2.7.14
cd /home/eld72413/DelMut/BAD_Mutations/Supporting
python VeP_to_Subs.py $VEP_OUTPUT_GZIP $OUTPUTFILE $OUTPUTDIR

```

#### Generate FASTA query files
Will start by testing 1 substitution region first

I used the gffread utility
```bash
module load gffread/0.9.12-foss-2016b

GFF3=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3
FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
OUTPUTDIR=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/FASTA_test

#Started with 1 representative sequence-
awk '{if ($3=="CDS") {print $0}}' Ha412HOv2.0-20181130.gff3 | grep mRNA:Ha412HOChr17g0858291

# test gffread on same region
gffread $GFF3 -g $FASTA -r Ha412HOChr17:205476753..205477742 -x ${OUTPUTDIR}/Ha412HOChr17g0858291_TEST.fasta

#count nucleotides
grep -v "^>" Ha412HOChr17g0858291_TEST.fasta | grep -Eo '[[:alnum:]]' | wc -l #990
```

gffread info:
`gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] 
 [-o <outfile.gff>] [-t <tname>] [-r [[<strand>]<chr>:]<start>..<end> [-R]]
 [-CTVNJMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]
 [-i <maxintron>] `

-s  #<seq_info.fsize> is a tab-delimited file providing this info
    # for each of the mapped sequences:
    # <seq-name> <seq-length> <seq-description>
    # (useful for -A option with mRNA/EST/protein mappings)
-C # coding only: discard mRNAs that have no CDS feature
-x # write a fasta file with spliced CDS for each GFF transcript
-W  #	for -w and -x options, also write for each fasta record the exon coordinates projected onto the spliced sequence
 -y   # write a protein fasta file with the translation of CDS for each record

this was the first way I tried to do it:
```bash
GFF3_file=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3
region=mRNA:Ha412HOChr17g0858291

Interval=$(awk '{if ($3=="CDS") {print $0}}' $GFF3_file | grep $region | awk '{print $1":"$4"-"$5}')

REF_FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
OUTPUTDIR=/scratch/eld72413/NSFproj/VEP/NewOutputOct2020/FASTA_test

module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

gatk FastaReferenceMaker \
	-R "$REF_FASTA" \
	-O "$OUTPUTDIR/Ha412HOChr17g0858291.fasta" \
	-L "${Interval}"
# I need to change the header on the FASTA to match...?
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


# location of genome alignments on MSI:
/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps/BAD_Mutations_Genomes

