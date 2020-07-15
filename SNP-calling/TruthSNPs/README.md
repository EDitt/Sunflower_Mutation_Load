# 10k SNPs

This directory contains the process used to map 10k SNPs from Bachlava et al. (2012): https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029814#s4 to new HA412v.2 genome  
Credit to Chaochih Liu for help with this process: (see https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/50k_9k_BOPA_SNP/README.md)

### Navigation: Jump to Section

---

## Data

A 10k SNP Illumina Infinium SNP array was created by aligning short-read transcriptome data against a reference transcriptome made from assembled long-read ESTs from several elite sunflower lines. The top 10,640 SNPs were selected based on coverage, minimum allele frequency, mapping quality and annotations (putative introns or non-genic regions discarded).  
For details see Bachlava et al. (2012) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0029814#s4  

This SNP array was used to genotype four mapping populations to construct a consensus linkage map: see Bowers et al. (2012) https://pubmed.ncbi.nlm.nih.gov/22870395/  

SNPs were also genotyped on a diverse collection of 271 sunflower lines (mostly overlapping with lines used to call SNPs in the current experiment). The total number of readily scorale, bi-allelic SNPs was 5,788, out of which 5,359 had MAF >/= 10%  
See Mandel et al. (2013): https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003378

---

## Data Exploration

Contextual sequences for the SNPs identified in Bachlava et al. (2012) were taken from supplementary file s006 and saved as a tab-delimited .txt file  
This file was converted to FASTA format using the following code:

```bash
awk 'NR >1 { print ">"$2"\n"$5 }' Bachlava_subset.txt > ContextualSeqs.fasta
```

I also replaced the polymorphism syntax (e.g. "[A/C]") with ambiguous characters designated by IUPAC codes
A/C - 1134; T/C - 4322; A/G - 3998; T/G - 1186 (No A/T or C/G polymorphisms, as these were removed)
```bash
sed -i 's|\[A/C]|M|g' ContextualSeqs.fasta
sed -i 's|\[T/C]|Y|g' ContextualSeqs.fasta
sed -i 's|\[A/G]|R|g' ContextualSeqs.fasta
sed -i 's|\[T/G]|K|g' ContextualSeqs.fasta
```

Index reference  
Adopted from Chaochih's script: https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/prep_reference/make_index_pseudo_bowtie2.sh  

```bash
qsub RefBuild.sh
```

Alignment
Adopted from Chaochih's script: https://github.com/MorrellLAB/morex_reference/blob/master/morex_v2/prep_reference/check_by_aligning_bowtie2_BOPA.sh

```bash
qsub
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
grep "^>" ContextualSeqs_subset.fasta | wc -l #5360
```

I also replaced the polymorphism syntax using the same code as above  
After subsetting: A/C - 507; T/C - 2208; A/G - 2115; T/G - 529  

I then ran the alignment again with this new set.



---

## Data Preparation

#### BLAST database

First, donwloaded taxonomy information:
```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
# Extract files
tar -xzvf taxdb.tar.gz
```

Created blast database:
```bash
module load BLAST+/2.10.0
makeblastdb -in Ha412HOv2.0-20181130.fasta -dbtype nucl -parse_seqids -title "Ha412HOv2_DB"
#check
blastdbcheck -db Ha412HOv2.0-20181130.fasta
```

blast search for the SNPs against the reference database
```bash
blastn -db Ha412HOv2.0-20181130.fasta -query /scratch/eld72413/SNParray/ContextualSeqs_subset.fasta -out /scratch/eld72413/SNParray/Blast_Mandel_subsetSNPs.out
```

#### Create Illumina Lookup Table  
This is a two-column, headerless table that has a SNP ID and contextual sequence in with the SNP in brackets and is a required input for SNP-utils  
I again subsetted so I only used the SNPs genotyped in Mandel et al. (2013)

```bash
for i in `awk -F"\t" 'NR >1 { print $1 }' Mandel_TableS2.txt `; do
	awk -v OFS='\t' -v var="$i" 'NR >1 {if ($2 == var) { print $2, $5 }}' Bachlava_subset.txt >> LookupTable_MandelSub.txt
done
```

#### Format Genetic Map
The creation of the genetic map is described in Bowers et al. (2012). I saved the first sheet "combined map" from supplementary file S1 from this paper as a tab-delimited .txt file
```bash
```