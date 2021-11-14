# 10k SNPs

This directory contains the process used to map 10k SNPs from Bachlava et al. (2012): https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029814#s4 to new HA412v.2 genome  
Credit to Chaochih Liu for help with this process: (see https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/README.md)

### Navigation: Jump to Section

- [Data](#data)
- [Data Preparation](#data-preparation)
- [Data Exploration](#data-exploration)
  - [Alignment Results](#bowtie2-alignment-results)
- [SNP-Utils](#snp-utils)
- [Results Summary](#Result_summary)
- [Clean-up](#clean-up)

---

## Data

A 10k SNP Illumina Infinium SNP array was created by aligning short-read transcriptome data against a reference transcriptome made from assembled long-read ESTs from several elite sunflower lines. The top 10,640 SNPs were selected based on coverage, minimum allele frequency, mapping quality and annotations (putative introns or non-genic regions discarded).  
For details see Bachlava et al. (2012) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029814#s4  

This SNP array was used to genotype four mapping populations to construct a consensus linkage map: see Bowers et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22870395/  
The consensus map contains 10,083 loci, and includes 1512 PCR-based loci (not from SNP array). 783 mapped to >1 loci (762 to 2 different locations and 21 to 3 different locations)

SNPs were also genotyped on a diverse collection of 271 sunflower lines (mostly overlapping with lines used to call SNPs in the current experiment). The total number of readily scorale, bi-allelic SNPs was 5,788, out of which 5,359 had MAF >/= 10%  
See Mandel et al. (2013): https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003378

---
## Data Preparation

#### Contextual Sequences
Contextual sequences for the SNPs identified in Bachlava et al. (2012) were taken from supplementary file s006 and saved as a tab-delimited .txt file  
This file was converted to FASTA format using the following code:

```bash
awk 'NR >1 { print ">"$2"\n"$5 }' Bachlava_subset.txt > ContextualSeqs.fasta

# re-created using only the SNPs mapped uniquely and used in the Bowers et al. (2012) genetic map (see below)
awk 'FNR==NR{arr[$0];next}($2 in arr)' Unique_Assay_Names.txt Bachlava_subset.txt |  awk '{ print ">"$2"\n"$5 }' > ContextualSeqs_UniqMap.fasta

```

I also replaced the polymorphism syntax (e.g. "[A/C]") with ambiguous characters designated by IUPAC codes  
Number of polymorphisms: A/C - 1134; T/C - 4322; A/G - 3998; T/G - 1186 (No A/T or C/G polymorphisms, as these were removed)
```bash
sed -i 's|\[A/C]|M|g' ContextualSeqs.fasta
sed -i 's|\[T/C]|Y|g' ContextualSeqs.fasta
sed -i 's|\[A/G]|R|g' ContextualSeqs.fasta
sed -i 's|\[T/G]|K|g' ContextualSeqs.fasta
```

#### Format Genetic Map
The creation of the genetic map is described in Bowers et al. (2012). I saved the first sheet "combined map" and the last sheet "assays that map to > 1 loci" from supplementary file S1 from this paper as separate tab-delimited .txt files.

Combined Map
```bash
# clean-up- remove top row and un-needed columns
awk 'NR >1 {print $0}' MapData_Bowers2012_FileS1.txt  | cut -f10-13 > MapData_Bowers2012_reduced.txt
# check numbers:
wc -l MapData_Bowers2012_reduced.txt #10084 (10,083 markers)
awk -F'\t' '$2!=""' MapData_Bowers2012_reduced.txt | wc -l #8572 (8571 without the 1512 PCR-based markers)
```

Assays that map to >1 loci
```bash
awk 'NR >2 {print $0}' MapData_Bowers2012_FileS1_multimappers.txt  | cut -f11  | sort -u > SNP_MultiMappers.txt
wc -l SNP_MultiMappers.txt #783

# how many markers are *not* multi-mappers?
grep -wvf SNP_MultiMappers.txt MapData_Bowers2012_reduced.txt | wc -l #8497 (including header line)
```

```bash
# genetic map for all SNPs that map uniquely to genome

awk 'FNR==NR{arr[$0];next}!($2 in arr)' SNP_MultiMappers.txt MapData_Bowers2012_reduced.txt | awk -F'\t' '$2!="" && NR >1 {print $3,$1,$4,"-"FNR}' OFS='\t' > SNP_Genetic_Map_Unique.txt

# there was one extra zero in the names from this paper compared to the names in the other two
sed -i 's|SFW0|SFW|g' SNP_Genetic_Map_Unique.txt

# check
wc -l SNP_Genetic_Map_Unique.txt #6984 (8496 unique mappers minus 1512 PCR-based markers)

# saved this list to subset others
awk '{print $2}' SNP_Genetic_Map_Unique.txt > Unique_Assay_Names.txt

# how many of these were genotyped in Mandel et al. (2013)?
awk 'FNR==NR{arr[$0];next}($1 in arr)' Unique_Assay_Names.txt Mandel_TableS2.txt | wc -l #5359 (all that were genotyped were in this unique set)

```

#### Create Illumina Lookup Table  
This is a two-column, headerless table that has a SNP ID and contextual sequence in with the SNP in brackets and is a required input for SNP-utils  

```bash
# for all SNPs
awk -v OFS='\t' 'NR >1 { print $2, $5 }' Bachlava_subset.txt > LookupTable_All.txt

# for uniquely mapped SNPs
awk 'FNR==NR{arr[$0];next}($2 in arr)' Unique_Assay_Names.txt Bachlava_subset.txt | awk -v OFS='\t' '{ print $2, $5 }' > LookupTable_MapUniq.txt

# for genotyped subset
for i in `awk -F"\t" 'NR >1 { print $1 }' Mandel_TableS2.txt `; do
  awk -v OFS='\t' -v var="$i" 'NR >1 {if ($2 == var) { print $2, $5 }}' Bachlava_subset.txt >> LookupTable_MandelSub.txt
done
```


#### Summary:
Bachlava et al. (2012) designed an array that included 10,640 SNPs (though only 9,480 were included due to 10.9% manufacturing loss)

Bowers et al. (2012) mapped 10,083 loci, including 1512 PCR-based markers (not from array). Of the 8571 array-based markers used, 6984 mapped to only 1 location in the genome.

Mandel et al. (2012) genotyped 5,359 array markers in 271 lines, All of these were in the "uniquely mapped" set from Bowers et al. (2012)

---

## Data Exploration

#### Mapping
Mapped contextual sequences to the new reference:
1.) Index reference  
Adopted from Chaochih's script: https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/prep_reference/make_index_pseudo_bowtie2.sh  

```bash
qsub RefBuild.sh
```

2.) Alignment
Adopted from Chaochih's script: https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/prep_reference/check_by_aligning_bowtie2_BOPA.sh

```bash
qsub AlignBowtie.sh
```

I will also subset these 10,640 SNPs to only include the 5,359 that were genotyped in Mandel et al. (2013)  

Supplementary Table 2 (associated mapping results for branching) from Mandel et al. (2013) was saved as a tab-delimited .txt file  
I subsetted the 10,640 list from Bachlava to the 5,359 based on the SNP names and made a new fasta file with this subset
```bash
for i in `awk -F"\t" 'NR >1 { print $1 }' Mandel_TableS2.txt `; do
	awk -v OFS='\t' -v var="$i" 'NR >1 {if ($2 == var) { print ">"$2"\n"$5 }}' Bachlava_subset.txt >> ContextualSeqs_subset.fasta
done
```
```bash
#check
grep "^>" ContextualSeqs.fasta | wc -l #10640
grep "^>" ContextualSeqs_UniqMap.fasta | wc -l #6984
grep "^>" ContextualSeqs_subset.fasta | wc -l #5360
```

I also replaced the polymorphism syntax using the same code as above  
```bash
# In dir: /scratch/eld72413/SAM_seq/Recombination
sed -i 's|\[A/C]|M|g' ContextualSeqs_UniqMap.fasta
sed -i 's|\[T/C]|Y|g' ContextualSeqs_UniqMap.fasta
sed -i 's|\[A/G]|R|g' ContextualSeqs_UniqMap.fasta
sed -i 's|\[T/G]|K|g' ContextualSeqs_UniqMap.fasta
```
Number of polymorphisms after subsetting to only include genotyped results: A/C - 507; T/C - 2208; A/G - 2115; T/G - 529  

I then ran the alignment again with this new set.

### Bowtie2 Alignment Results
All sequences; N= 10640
```bash
10640 reads; of these:
  10640 (100.00%) were unpaired; of these:
    1672 (15.71%) aligned 0 times
    6703 (63.00%) aligned exactly 1 time
    2265 (21.29%) aligned >1 times
84.29% overall alignment rate
```
Subset that mapped uniquely in Bowers et al. (2012); N= 6984
```bash
6984 reads; of these:
  6984 (100.00%) were unpaired; of these:
    910 (13.03%) aligned 0 times
    4741 (67.88%) aligned exactly 1 time
    1333 (19.09%) aligned >1 times
86.97% overall alignment rate
```

Subset that were genotyped in Mandel et al. (2013); N= 5359
```bash
5359 reads; of these:
  5359 (100.00%) were unpaired; of these:
    564 (10.52%) aligned 0 times
    4003 (74.70%) aligned exactly 1 time
    792 (14.78%) aligned >1 times
89.48% overall alignment rate
```

Summary stats:

```bash
module load SAMtools/1.10-GCC-8.3.0
samtools flagstat Bachlava_Contextual_HA412v2_bowtie2.sam
10640 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
8968 + 0 mapped (84.29% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
Summary stats for subset that mapped uniquely
```bash
6984 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
6074 + 0 mapped (86.97% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Get summary stats for subset that were genotyped
```bash
samtools flagstat Contextual_Subset_HA412v2_bowtie2.sam 
5359 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
4795 + 0 mapped (89.48% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

---

## SNP-Utils 

#### Create BLAST database

First, downloaded taxonomy information:
```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# Extract files
tar -xzvf taxdb.tar.gz
```

Created blast database:
```bash
#module load BLAST+/2.10.0
module load BLAST+/2.11.0-gompi-2020b
# In dir: /scratch/eld72413/SunflowerGenome
makeblastdb -in Ha412HOv2.0-20181130.fasta -dbtype nucl -parse_seqids -title "Ha412HOv2_DB"

# -taxid 4232 do I need this flag?
#check
blastdbcheck -db Ha412HOv2.0-20181130.fasta
```
(initially got an error when I didn't extract the taxonomy info from the BLAST database in the same directory as the reference genome)

blast search for the SNPs against the reference database
```bash
blastn -db Ha412HOv2.0-20181130.fasta -query /scratch/eld72413/SNParray/ContextualSeqs_UniqMap.fasta -out /scratch/eld72413/SNParray/Blast_UniqueSNPs.out

# redo (diff directory)
blastn -db Ha412HOv2.0-20181130.fasta -query /scratch/eld72413/SAM_seq/Recombination/ContextualSeqs_UniqMap.fasta -out /scratch/eld72413/SAM_seq/Recombination/Blast_UniqueSNPs.out
# lots of this message: FASTA-Reader: Ignoring invalid residues at position(s): On line 13968: 37, 39, 41
```

#### Run SNP-Utils
First, install dependencies  
Listed here: https://github.com/mojaveazure/SNP_Utils#dependencies
```bash
#module load BLAST+/2.10.0
# redo
module load BLAST+/2.11.0-gompi-2020b
module load Anaconda3/2020.02

# use PyPi to install biopython, Beautiful Soup 4, Overload, Ixml
pip install PyPi
pip install biopython
pip install beautifulsoup4
pip install overload
pip install lxml

# to update pip: /apps/eb/Anaconda3/2020.02/bin/python -m pip install --upgrade pip
```

Step 1: Config subroutine (configure BLAST search)
```bash
#GENOME=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
GENOME=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta
cd /home/eld72413/DelMut/SNP_Utils/
./snp_utils.py CONFIG -d ${GENOME} -k -i 90 -c /scratch/eld72413/SNParray/SNPutils/blast_MapUniqueSNP_idt90
```

```bash
Validating config
Searching for proper database files for Ha412HOv2.0-20181130.fasta in /scratch/eld72413/Ha412HOv2.0
Setting option database with value /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta
Setting option evalue with value 0.1
Setting option max_hits with value 3
Setting option max_hsps with value 3
Setting option identity with value 90.0
Setting option keep_query with value True

Config file can be found at /scratch/eld72413/SNParray/SNPutils/blast_MapUniqueSNP_idt90
```

Step 2: Run SNP-Utils BLAST
```bash
# define variables
#LOOKUP_TABLE=/scratch/eld72413/SNParray/LookupTable_MapUniq.txt
LOOKUP_TABLE=/scratch/eld72413/SAM_seq/Recombination/LookupTable_MapUniq.txt
#GENETIC_MAP=/scratch/eld72413/SNParray/SNP_Genetic_Map_Unique.txt
GENETIC_MAP=/scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
OUT_PREFIX=/scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90

cd /home/eld72413/DelMut/SNP_Utils
./snp_utils.py BLAST -l ${LOOKUP_TABLE} -c /scratch/eld72413/SNParray/SNPutils/blast_MapUniqueSNP_idt90 -b -m ${GENETIC_MAP} -d -t 100000 -o ${OUT_PREFIX}
```

```bash
Filtering SNPs by hit chromsome/contig
Using genetic map /scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
Filtering SNPs with a minimum distance threshold of 100000
Filtering SNPs by relative location on the genetic map
Using genetic map /scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
Parsing blast database /scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta
Found 25090 chromosomes
Filtering 6539 SNP IDs
Writing 6539 SNPs to /scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90.vcf
Removing masked SNPs that were actually found
Writing 3 masked SNPs to /scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90_masked.vcf
Writing 134 failed SNPs to /scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90_failed.log
```

### Result summary
```bash
cd /scratch/eld72413/SNParray/SNPutils/
# Total SNPs
grep -v "#" MapUniqueSNP_idt90.vcf | cut -f 3 | wc -l # 6539
# Total unique SNPs
grep -v "#" MapUniqueSNP_idt90.vcf | cut -f 3 | sort -u | wc -l #6539
# Total number of duplicates
grep -v "#" MapUniqueSNP_idt90.vcf | cut -f 3 | sort | uniq -c | sort -n -r | grep -vw "1" | wc -l #0

```
Running `snp_utils.py` for 90idt, I got:
- 6,539 unique SNPs
- 0 SNPs with duplicates
- 3 masked SNPs
- 134 failed SNPs

### Clean up

I need to do two things-   
1.) remove SNPs that are located in contigs less than 10kbp. I did not call SNPs in these regions.   
2.) Convert length positions to chromosome names. For some reason SNP-utils listed the chromosome lengths instead of names in the .vcf.

I need to do both of these together because the lengths of chromosomes are repeated in the small contigs, so they are unable to be distinguished

How many SNPs were in small contigs (smaller than 10kbp)?
```bash
VCF="/scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90.vcf"

grep -v "#" $VCF | awk '{print $1}' | cut -c 5- | awk '{if ($1 < 10000) {print $0}}' | wc -l # Answer: 15 (these will be removed)

```
  
Use reference dictionary to find and replace:  
```bash
# First make a list of chromosome names and lengths:
#DICT="/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.dict"
DICT="/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.dict"

# make a list of chromosome names + lengths
awk -F "[\t,:]" 'NR > 1 {$1=$1; print $3,$5}' $DICT > Chrom_Names_Len.txt

# select only contigs greater than 10,000
awk '{if ($2 > 10000) {print $0}}' Chrom_Names_Len.txt > Chrom_Names_Len_Over10kbp.txt
wc -l Chrom_Names_Len_Over10kbp.txt #590

```

Subset to get only those represented in VCF (5 lengths are represented more than once in this set)
```bash
VCF="/scratch/eld72413/SNParray/SNPutils/MapUniqueSNP_idt90.vcf"

grep -v "#" $VCF | awk '{print $1}' | sort -u | wc -l #36

# lengths represented in VCF
VCFLengths=$(grep -v "#" $VCF | awk '{print $1}' | sort -u )

# make a new list that only contains the chromosomes represented in the VCF
# this is the union of a.) in VCF and b.) over 10kbp
for i in $VCFLengths; do
  awk -v var="${i#len=}" '{if ($2 == var) {print $1,"len="$2}}' Chrom_Names_Len_Over10kbp.txt >> Chrom_Names_Len_InVCF_Over10kbp.txt
done

wc -l Chrom_Names_Len_InVCF_Over10kbp.txt #21 (17 chromosomes + 4 contigs over 10kbp) -> 36 minus the 15 small contigs
```

Make a new VCF file with contigs renamed and excluding small contigs
```bash
cp $VCF ./MapUniqueSNP_idt90_rename.vcf

while read line; do
  ChromName=$(echo $line | cut -d " " -f1)
  ChromLen=$(echo $line | cut -d " " -f2)
  sed -i 's|'^"${ChromLen}"'\b|'"${ChromName}"'|g' MapUniqueSNP_idt90_rename.vcf
done < Chrom_Names_Len_InVCF_Over10kbp.txt

#check
grep "^len=" MapUniqueSNP_idt90_rename.vcf | wc -l #15

```

Remove 15 small contigs that could not be re-named
```bash
grep -v "^len=" MapUniqueSNP_idt90_rename.vcf > MapUniqueSNP_idt90_rename_rmContigs.vcf 

# how many markers left?
grep -v "#" MapUniqueSNP_idt90_rename_rmContigs.vcf | wc -l #6524
```

Sort by Position
```bash
qsub -I -q s_interq -l walltime=8:00:00 -l nodes=1:ppn=8 -l mem=22gb

module load picard/2.16.0-Java-1.8.0_144
PICARD_JAR=/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
GATK_JAR=/usr/local/apps/eb/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk

TMP=/scratch/eld72413/SAM_seq/results2/Temp
DIR=/scratch/eld72413/SNParray

gatk SortVcf \
--TMP_DIR ${TMP} \
-I ${DIR}/MapUniqueSNP_idt90_rename_rmContigs.vcf \
-SD /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.dict \
-O ${DIR}/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf
```


gzip for QC with bcftools
```bash
bgzip $Truth_Set
tabix MapUniqueSNP_idt90_rename_rmContigs.vcf.gz
```