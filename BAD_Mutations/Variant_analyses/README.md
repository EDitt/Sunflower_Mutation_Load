# dSNP Analysis

## Navigation: Jump to Section

- [Annotation Classes](#parse-vep-table-and-predictions-output)
- [Polarize Ancestral State](#polarize-ancestral-state)

- [Subset VCF](#subset-vcf-by-consequence)
- [Subset VCF](#subset-vcf-by-consequence)
- [Site Frequency Spectra](#site-frequency-spectra)
- [Genomic Patterns](#genomic-patterns)
	- [Recombination across genome](#recombination-across-genome)
	- [Binning SNP classes](#binning-snp-classes)
- 

---

## Parse VeP Table and Predictions output
See: `Variant_class_numbers.md` for commands used to parse the prediction output and VeP output to get the variants for each class

## Polarize Ancestral State

Created ancestral fasta file with ANGSD. See: ANGSD directory

See `Polarize_SNPs.md` for process of creating ancestral state table.
Output a table with ancestral state calls at 13,524,630 variant positions: `/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt`


## Site Frequency Spectra

First, parsed VeP table and dSNP predictions to get the positions of SNPs in each frequency class (see `Variant_class_numbers.md`)

Next, get allele frequency information for all alleles
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.13-GCC-8.3.0
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' ${vcf} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt

# take out the sites with allele count of 1 (heterozygous in only 1 sample due to sample pooling)
Pos_info=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt
wc -l  $Pos_info # 37,129,915
awk '{if ($7<0.05) {print $0}}' $Pos_info | wc -l # 28,667,994 in this frequency class

awk '{if ($5==1) {print $0}}' $Pos_info | wc -l # 7,805,470
awk '{if ($5==2) {print $0}}' $Pos_info | wc -l # 7,625,046
awk '{if ($5>1) {print $0}}' $Pos_info > All_alleleFreqInfo_noSingleton.txt # 29,324,445

# what about cases where reference exhibits the minor allele only in heterozygote form?

# use R to get histogram bin information
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses
module load R/4.0.0-foss-2019b
Rscript "SFS_Info.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo_noSingleton.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/MAF_Bins.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedFreq_Bins.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/FrequencyInfo.RData" # this file will be used for other purposes

```

However, to look at derived dSNPs, I need to separate the deleterious alleles in the reference and alternate states (i.e. if a dSNP is in the alternate state, but the reference allele is derived, this allele would not be counted as derived deleterious)

Proportions in each frequency class:
```bash
MAF=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/MAF_Bins.txt 
Derived_freq=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedFreq_Bins.txt
Derived_dSNPfreq=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Derived_dSNP_freqbins.txt

awk '{if ($3=="AllDel") {print $0}}' $MAF
awk '{if ($3=="NonCoding") {print $0}}' $MAF
awk '{if ($3=="StopLostGained") {print $0}}' $MAF
awk '{if ($3=="Synonymous") {print $0}}' $MAF
awk '{if ($3=="Tolerated") {print $0}}' $MAF

awk '{if ($3=="NonCoding") {print $0}}' $Derived_freq
awk '{if ($3=="StopLostGained") {print $0}}' $Derived_freq
awk '{if ($3=="Synonymous") {print $0}}' $Derived_freq
awk '{if ($3=="Tolerated") {print $0}}' $Derived_freq

```

#### Deleterious SNPs
Compiled Report was used to identify variant classes (see `2.BAD_Mutations/Post_processing.md` and positions files created in `Variant_class_numbers.md`)

First, subset vcf file to make separate vcfs of deleterious (for both reference and alternate alleles)
```bash
cd /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations

# deleterious in reference
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Reference_DelPositions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Refdeleterious' Subset_vcf.sh # Submitted batch job 4232837

grep -v "#" SAM_Refdeleterious.vcf | wc -l # 11796

# deleterious in alternate
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Alternate_DelPositionsNoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_Altdeleterious' Subset_vcf.sh # Submitted batch job 4232853

grep -v "#" SAM_Altdeleterious.vcf | wc -l # 76016

# use bcftools to get stats for each allele
module load BCFtools/1.13-GCC-8.3.0

# need to bgzip vcf files:
module load HTSlib/1.10.2-GCC-8.3.0
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf.gz

# derived deleterious in alternate
bgzip -c --threads 4 /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz
tabix -p vcf /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf.gz

# reference deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Refdeleterious.vcf.gz > AlleleFreqs/SAM_Refdeleterious_AC_AN.txt
# alternate deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Altdeleterious.vcf.gz > AlleleFreqs/SAM_Altdeleterious_AC_AN.txt

# N_ALT : number of alternate alleles
# AC : count of alternate alleles
# AN : number of alleles in called genotypes

### need to calculate frequency of deleterious allele (for deleterious reference it's AN - AC/AN; for deleterious alternate it's AC/AN)
```

Use R to wrangle
```R
#module load R/4.0.0-foss-2019b
#R

dSNPNums <- function (Alt_file, Ref_file, AncestralState_file) {
	Alt_dSNP <- read.table(Alt_file, sep = "\t", header=FALSE)
	colnames(Alt_dSNP) <- c("Chromosome", "Position", "Ref_allele", "Del_allele",
		"Num_Alt_alleles", "Num_alleles")
	Alt_dSNP$dSNP_freq <- Alt_dSNP$Num_Alt_alleles / Alt_dSNP$Num_alleles
	Ref_dSNP <- read.table(Ref_file, sep = "\t", header=FALSE)
	colnames(Ref_dSNP) <- c("Chromosome", "Position", "Del_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles")
	Ref_dSNP$dSNP_freq <- (Ref_dSNP$Num_alleles - Ref_dSNP$Num_Alt_alleles) / Ref_dSNP$Num_alleles
	All_dSNP_freq <- rbind(Ref_dSNP[,c(1:3,7)], Alt_dSNP[,c(1,2,4,7)])
	All_dSNP_freq$Chromosome <- as.factor(All_dSNP_freq$Chromosome)
	All_dSNP_freq <- with(All_dSNP_freq, All_dSNP_freq[order(Chromosome, Position),])
	Polarized <- read.table(AncestralState_file)
	colnames(Polarized) <- c("Chromosome", "Position", "AncestralAllele")
	dSNP_summary_AA <- merge(All_dSNP_freq, Polarized, by=c("Chromosome", "Position"), all.x=TRUE)
	dSNP_summary_AA$Cat <- ifelse(dSNP_summary_AA$AncestralAllele==dSNP_summary_AA$Del_allele, "Ancestral_dSNP",
	ifelse(dSNP_summary_AA$AncestralAllele!=dSNP_summary_AA$Del_allele,"Derived_dSNP", NA))
  	return(dSNP_summary_AA)
}

dSNP_summary <- dSNPNums("SAM_Altdeleterious_AC_AN.txt", "SAM_Refdeleterious_AC_AN.txt", "/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt")

#duplicates?
length(unique(paste0(dSNP_summary$Chromosome,"_", dSNP_summary$Position))) # 87794 (out of 87812) 18 duplicated
# how many of these are polarized?
aggregate(dSNP_summary$Position, by=list(dSNP_summary$Cat), length)
# ancestral dSNP= 6141; derived dSNP= 56360 (total 62,501) <- still have 14 duplicate positions

derived_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
Derived_dSNP_histogram <- hist(dSNP_summary[which(dSNP_summary$Cat=="Derived_dSNP"), "dSNP_freq"],
	plot = FALSE, breaks=derived_breaks)
sum(Derived_dSNP_histogram$counts) # 56360

# use "Hist_bins" function from SFS_Info.R
Derived_dSNP <- Hist_bins(dSNP_summary[which(dSNP_summary$Cat=="Derived_dSNP"),],
	derived_breaks, "dSNP_freq", "dSNP")
colnames(Derived_dSNP)[1] <- "breaks"
write.table(Derived_dSNP, file="/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Derived_dSNP_freqbins.txt", sep = "\t", quote=FALSE, row.names=FALSE)

```


----- 

Continued with `dSNP_DerAncBINS.R` <- for binning across genome

----- 



## Genomic Patterns

### Recombination across genome

Using John Bowers genetic map information, I used the file I manipulated for the truth SNPs: `SNP_Genetic_Map_Unique.txt` for the cM distance of the markers. This file contains 6984 markers that mapped uniquely. 
- The first column is the linkage group (chromosome number), second column is the locus name, third column is the distance in cM.

I also used the vcf file that I created with SNPutils- `MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf` where I mapped the SNPs to the new genome build. (N=6523)
- Here the first column is the chromosome, 2nd is the position in bp, third is Locus name

##### Create Recombination Data Files
```bash
GeneticMap=/scratch/eld72413/SAM_seq/Recombination/SNP_Genetic_Map_Unique.txt
RemappedVCF=/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf

cd /scratch/eld72413/SAM_seq/Recombination

grep -v "#" $RemappedVCF | awk '{print $1,"\t",$2,"\t",$3}' > SNParray_BPpositions.txt

awk '{print $1,"\t",$2,"\t",$3}' $GeneticMap > SNParray_cMpositions.txt


# I need to change the chromosome names to match:

awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' # need to remove the four contigs
awk '{print $1}' SNParray_BPpositions.txt | sort -u | awk 'NR > 4 {print $0}' | cat -n > ChromNames.txt 


while read line; do
	OldChromName=$(echo $line | awk '{print $1}')
	NewChromName=$(echo $line | awk '{print $2}')
	echo "$OldChromName to $NewChromName"
	sed -i 's/'^"${OldChromName}"'\b/'"${NewChromName}"'/g' SNParray_cMpositions.txt
done < ChromNames.txt
```

To graph Recombination patterns across genome, I used R (see Genomic_Patterns/Recombination.R script).

## Binning SNP classes across Genome


### Number of codons per window
```bash
GFF3=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.gff3
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

awk '{if ($3=="mRNA") {print $1,$4,$5}}' $GFF3 | wc -l # 72995

# use bedtools to find # bp (and # counts) for each mRNA feature in 10Mbp window- divide by 3 to get number of codons
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b $GFF3 -wo \
| awk '{if ($6=="mRNA") {print $0}}' \
| bedtools groupby -g 1,2,3,6 -c 13 -o sum,count \
| awk '{print $1,$2,$3,$4, $5,$6, $5/3}' > ${OUTPUTDIR}/mRNA_bpNumCodonCounts.txt

```

### Number of different SNP classes per window
Using R object saved from SFS_Info.R script
```R
# tmux new -s SNP_bins
# srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# module load R/4.0.0-foss-2019b
load("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/FrequencyInfo.RData") # object is FrequencyInfo

# function to bin number of each SNP class across genome
SNP_bins <- function(dataset, BinSize_bp, SNP_class_Name) {
	dataset <- dataset[order(dataset$Position),]
	Num_windows <- ceiling(max(dataset$Position)/BinSize_bp)
	dataset$bin <- cut(dataset$Position,seq(0,Num_windows*BinSize_bp,BinSize_bp))
	counts <- aggregate(dataset$Position,
                      by=list(dataset$bin), length, drop=FALSE)
	colnames(counts) <- c("Bin", paste0("Number_", SNP_class_Name))
	counts$BinList <- strsplit(as.character(counts$Bin), ",")
	counts$StartPos <- as.numeric(gsub('[(]', '', sapply(counts$BinList, "[", 1)))
	counts$EndPos <- as.numeric(gsub('[]]', '', sapply(counts$BinList, "[", 2)))
	return(counts[,c(4,5,2,1)])
}

# function to apply SNP_bins across all chromosomes
apply_SNP_bins <- function(dataset, BinSize_bp, SNP_class_Name) {
	dataset$Chromosome <- factor(dataset$Chromosome)
	dataset_chrom <- split(dataset, dataset$Chromosome)
	binCounts <- lapply(dataset_chrom, function(x) {
		SNP_bins(x, BinSize_bp, SNP_class_Name)
	})
	binCounts <- lapply(names(binCounts), function(x) {
		binCounts[[x]]["Chromosome"] <- x; return(binCounts[[x]])
	})
	binCounts_df <- do.call("rbind", binCounts)
	return(binCounts_df)
}

### first subset to only derived SNPs
FrequencyInfo_Derived <- lapply(FrequencyInfo, function(x) {
	subset(x, !is.na(x$Derived_Freq))
	})
names(FrequencyInfo_Derived)
names(FrequencyInfo_Derived)[2:5] # not using the deleterious SNPs from this list

Derived_Bins <- lapply(names(FrequencyInfo_Derived)[2:5], function(x) {
	apply_SNP_bins(FrequencyInfo_Derived[[x]], 10000000, x)
	})
names(Derived_Bins) <- names(FrequencyInfo_Derived)[2:5]

## for dSNPs:
setwd("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/AlleleFreqs")
# use dSNPNums function above
dSNP_summary <- dSNPNums("SAM_Altdeleterious_AC_AN.txt", "SAM_Refdeleterious_AC_AN.txt", "/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt")
dSNP_derived <- subset(dSNP_summary, Cat=="Derived_dSNP")

# (take out scaffold seq)
dSNP_Derived_Bins <- apply_SNP_bins(dSNP_derived[-grep("Ha412HOChr00c", dSNP_derived$Chromosome),], 
	10000000, "dSNP")

### merge all classes
All_info1 <- Reduce(function(x, y) merge(x, y, 
	by=c("StartPos", "EndPos", "Chromosome", "Bin"), all=TRUE), Derived_Bins)

### merge (without scaffold, N=325) with dSNP data
All_info <- merge(dSNP_Derived_Bins, 
	All_info1[-grep("Ha412HOChr00c", All_info1$Chromosome),], 
	by=c("StartPos", "EndPos", "Chromosome", "Bin"), all=TRUE)

All_info$StartPos <- as.integer(All_info$StartPos)
All_info$EndPos <- as.integer(All_info$EndPos)

### combine with number of codon info
NumCodons <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/mRNA_bpNumCodonCounts.txt",
	header=FALSE)
colnames(NumCodons) <- c("Chromosome", "StartPos", "EndPos", "feature", "Total_bases", "Total_features", "Num_codons")

All_info2 <- merge(All_info, NumCodons[,c(1:3,7)], by=c("StartPos", "Chromosome"))
All_info2 <- All_info2[order(All_info2$Chromosome, All_info2$StartPos),]

write.table(All_info2, file="/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/Derived_VariantNums.txt",
	row.names=FALSE, quote=FALSE, sep="\t")

```


#### Older way I binned across genome:
Make windows to bin count of each class of SNPs-
dSNPs, sSNPs, Tolerated SNPs
```bash
# make genome file
cd /scratch/eld72413/SunflowerGenome
awk -v OFS='\t' {'print $1,$2'} "Ha412HOv2.0-20181130.fasta.fai" | head -17 > "GenomeFile.txt"

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
mkdir GenomicBins

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BEDTools/2.30.0-GCC-8.3.0

vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins

### dSNPs per window
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Refdeleterious.vcf -c > ${out_dir}/RefDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_Altdeleterious.vcf -c > ${out_dir}/AltDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_tolerated.vcf.gz -c > ${out_dir}/Tolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_synonymous.vcf.gz -c > ${out_dir}/Synonymous_10MbCounts.txt

### do the same for the derived alleles
vcf_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs
out_dir=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/GenomicBins

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedDeleterious.vcf -c > ${out_dir}/RefDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedTolerated.vcf -c > ${out_dir}/RefDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_RefDerivedSynonymous.vcf -c > ${out_dir}/RefDerivedSynonymous_10MbCounts.txt


##
bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedDeleterious.vcf -c > ${out_dir}/AltDerivedDeleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedTolerated.vcf -c > ${out_dir}/AltDerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SunflowerGenome/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b ${vcf_dir}/SAM_AltDerivedSynonymous.vcf -c > ${out_dir}/AltDerivedSynonymous_10MbCounts.txt

```

### Derived SNP Patterns for Different Classes of Germplasm
Count Number of Alt/Ref Alleles per genotype
(using bcftools)
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_RefDerivedDeleterious.vcf |grep "PSC" > Stats/RefDerivedDeleterious_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedDeleterious.vcf | grep "PSC" > Stats/AltDerivedDeleterious_perSampleCounts.txt

bcftools stats -s - SAM_RefDerivedTolerated.vcf | grep "PSC" > Stats/RefDerivedTolerated_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedTolerated.vcf | grep "PSC" > Stats/AltDerivedTolerated_perSampleCounts.txt

bcftools stats -s - SAM_RefDerivedSynonymous.vcf | grep "PSC" > Stats/RefDerivedSynonymous_perSampleCounts.txt
bcftools stats -s - SAM_AltDerivedSynonymous.vcf | grep "PSC" > Stats/AltDerivedSynonymous_perSampleCounts.txt

### need total number of genotypes:
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
bcftools stats -s - $vcf | grep "PSC" > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllVariantStats.txt

### need total number of derived genotypes:

vcf=
bcftools stats -s - $vcf | grep "PSC" > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllDerived_VariantStats.txt

module load R/4.0.0-foss-2019b
Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Germplasm/Derived_Variant_Numbers.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/Stats" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AllVariantStats.txt" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/DerivedAlleles/AlleleClassVCFs/AlleleNums_Geno.txt"


```
