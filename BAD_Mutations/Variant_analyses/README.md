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


## Site Frequency Spectrum

### Folded 

First, parsed VeP table and dSNP predictions to get the positions of SNPs in each frequency class (see `Variant_class_numbers.md`)

Next, get allele frequency information for all alleles
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.13-GCC-8.3.0
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' ${vcf} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt

# use R to get histogram bin information
cd
module load R/4.0.0-foss-2019b
Rscript "SFS_Info.R" "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles" "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt" "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleFreqBins.txt"
```

Use R to wrangle
```R
#module load R/4.0.0-foss-2019b
#R

Position_lists <- ImportTxts("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles")
Positions_annotate <- do.call("rbind", Position_lists)
FrequencyInfo <- SNP_freq("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt", Positions_annotate)
MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
MAF_histList <- lapply(names(FrequencyInfo), function(x) {
	Hist_bins(FrequencyInfo[[x]], MAF_breaks, "MAF", x)
	})
MAF_histogram_data <- do.call("rbind", MAF_histList)
write.table(MAF_histogram_data, "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleFreqBinsTESTSCRIPT.txt", sep = "\t", quote=FALSE, row.names=FALSE)


```

----- 

## Subset VCF by Consequence
Compiled Report was used to identify variant classes (see `2.BAD_Mutations/Post_processing.md`)

Created VCFs of the variants in the different allele classes (see `Subset_VCFbyConsequence.md`)



## Site Frequency Spectra
All dSNPs
```bash
module load BCFtools/1.13-GCC-8.3.0

# reference deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Refdeleterious.vcf.gz > AlleleFreqs/SAM_Refdeleterious_AC_AN.txt
# alternate deleterious
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n' SAM_Altdeleterious.vcf.gz > AlleleFreqs/SAM_Altdeleterious_AC_AN.txt

# N_ALT : number of alternate alleles
# AC : count of alternate alleles
# AN : number of alleles in called genotypes

### need to calculate frequency of deleterious allele (for deleterious reference it's AN - AC/AN; for deleterious alternate it's AC/AN)
# because dSNPs are polarized to some degree, I cannot compare full set to tolerated or synonymous 
```

Use R to wrangle
```R
#module load R/4.0.0-foss-2019b
#R

dSNPNums <- function (Alt_file, Ref_file) {
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
  	return(All_dSNP_freq)
}

dSNP_summary <- dSNPNums("SAM_Altdeleterious_AC_AN.txt", "SAM_Refdeleterious_AC_AN.txt")

dSNP_histogram <- hist(dSNP_summary$dSNP_freq, plot = FALSE)
#save(dSNP_histogram, file = "All_dSNPFreqData.RData")

# how many of these are polarized?
Polarized <- read.table("/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt")
colnames(Polarized) <- c("Chromosome", "Position", "AncestralAllele")

dSNP_summary_AA <- merge(dSNP_summary, Polarized, by=c("Chromosome", "Position"), all.x=TRUE)

#duplicates?
length(unique(paste0(dSNP_summary_AA$Chromosome,"_", dSNP_summary_AA$Position))) # 87794 (out of 87812) 18 duplicated

dSNP_summary_AA$Cat <- ifelse(dSNP_summary_AA$AncestralAllele==dSNP_summary_AA$Del_allele, "Ancestral_dSNP",
	ifelse(dSNP_summary_AA$AncestralAllele!=dSNP_summary_AA$Del_allele,"Derived_dSNP", NA))

aggregate(dSNP_summary_AA$Position, by=list(dSNP_summary_AA$Cat), length)
# ancestral dSNP= 6141; derived dSNP= 56360

# frequencey of derived & ancestral:
Derived_dSNP_histogram <- hist(dSNP_summary_AA[which(dSNP_summary_AA$Cat=="Derived_dSNP"), "dSNP_freq"],
	plot = FALSE)
Ancestral_dSNP_histogram <- hist(dSNP_summary_AA[which(dSNP_summary_AA$Cat=="Ancestral_dSNP"), "dSNP_freq"],
	plot = FALSE)
save(dSNP_histogram, Derived_dSNP_histogram, Ancestral_dSNP_histogram, 
	file = "All_dSNPFreqData.RData")
write.table(dSNP_summary_AA, file="dSNP_freq_summary.txt", row.names=FALSE, quote=FALSE)

```
Continued with `dSNP_DerAncBINS.R`



All SNPs across the three classes of variants: deleterious, tolerated, synonymous
```bash

```




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

### Binning SNP classes across Genome
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


