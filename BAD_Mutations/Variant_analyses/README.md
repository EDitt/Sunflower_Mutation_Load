# dSNP Analysis

## Navigation: Jump to Section

- [Polarize Ancestral State](#polarize-ancestral-state)
- [SNP Annotation Classes](#snp-annotation-classes)
- [Create SNP Table](#create-snp-table)
- [Site Frequency Spectra](#site-frequency-spectra)
- [Germplasm Patterns](#germplasm-patterns)



- [Genomic Patterns](#genomic-patterns)
	- [Recombination across genome](#recombination-across-genome)
	- [Binning SNP classes](#binning-snp-classes)
- 

---

The compiled dSNP report was merged with VeP output (missense variants only) and processed using criteria detailed in `${REPO_DIR}/2.BAD_Mutations/Post_processing.md`

Resulting file: `/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table`

## Polarize Ancestral State

Created ancestral fasta file with ANGSD. See: ANGSD directory

See `${REPO_DIR}/Variant_analyses/Polarize_SNPs.md` for process of creating ancestral state table.
Output a table with ancestral state calls at 13,524,630 variant positions: `/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt`

Merged ancestral states with dsnp_data.table file to get: `/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_Polarized.table`

## SNP Annotation Classes
Parse VeP Table and Predictions output to obtain annotation classes for all SNPs

See: `${REPO_DIR}/Variant_analyses/Variant_class_numbers.md` for commands used to parse the prediction output and VeP output to get the variants for each class
Tab-delimited lists of variants positions for each annotation class found in: `/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles`

## Create SNP Table
Table with SNP Frequency, Ancestral State, Genotype Calls, Annotation, for all SNPs

Use bcftools to output a table with Allele Position, Reference and Alternate alleles, alternate allele count, total allele count, alternate allele frequency
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' ${VCF} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt

# remove N's from ancestral state table to save memory
awk '{if ($3 != "N") {print $0}}' ${ANCESTRAL_STATE} > /scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt

# Use R script to combine information
Rscript "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/Variant_Table.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt" \
"/scratch/eld72413/SAM_seq/Polarized/AncestralStateCalls.txt" \
"/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info.txt"
```

## Site Frequency Spectra
Get frequency bins for graphing (folded and unfolded)

To look at derived dSNPs, I need to remove variants for which the derived variant is not the deleterious variant
```bash
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

## Folded SFS:
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${SNP_INFO}" \
"0.5" \
"0.05" \
"MAF" \
"/scratch/eld72413/SAM_seq/dSNP_results/MAF_Bins.txt"

## UnFolded SFS:
#	=> First remove derived frequency information for deleterious SNPs for variants at which the derived variant != the deleterious variant
awk '{if (($39 == "Alternate_deleterious" && $7 == "Ref_derived") ||  ($39 == "Reference_deleterious" && $7 == "Alt_derived")) {print $0}}' ${DSNP_DATA} | wc -l # 6144
awk 'BEGIN{FS=OFS="\t"}; {if (($39 == "Alternate_deleterious" && $7 == "Ref_derived") ||  ($39 == "Reference_deleterious" && $7 == "Alt_derived")) {print $8,$9}}' ${DSNP_DATA} | sort -V > ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt

# remove these lines from dSNP_data from Unfolded SFS calculations:
grep -w -Ff ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt ${SNP_INFO} | wc -l # 6141 (3 sites not present in the SNP_INFO table)
grep -v -w -f ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt ${SNP_INFO} | sort -V > ${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt

## Unfolded SFS:
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SFS_Info.R" \
"${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt" \
"1.0" \
"0.05" \
"Derived_Freq" \
"/scratch/eld72413/SAM_seq/dSNP_results/DerivedFreq_Bins.txt"

rm ${OUT_DIR}/IntermediateFiles/ToRemove_NotDerived.txt
```

See Plots.R for SFS plot

Do I want to remove singletons?
```bash
# take out the sites with allele count of 1 (heterozygous in only 1 sample due to sample pooling)
Pos_info=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/All_alleleFreqInfo.txt
wc -l  $Pos_info # 37,129,915
awk '{if ($7<0.05) {print $0}}' $Pos_info | wc -l # 28,667,994 in this frequency class

awk '{if ($5==1) {print $0}}' $Pos_info | wc -l # 7,805,470
awk '{if ($5==2) {print $0}}' $Pos_info | wc -l # 7,625,046
awk '{if ($5>1) {print $0}}' $Pos_info > All_alleleFreqInfo_noSingleton.txt # 29,324,445

# what about cases where reference exhibits the minor allele only in heterozygote form?
```

----- 

## Germplasm patterns
Spreadsheet `scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/LineKeywINFO.csv` contains grouping information for each genotype to match with Variant ID's from VCF

See `${REPO_DIR}/2.BAD_Mutations/PCA.md` for information on the creation of a PCA for genotypes in the SAM pop

Plotted Fst across the genome. See:

## Heterotic Group Differentiaton

First, I will look at the proportion of shared versus private SNPs for deleterious and synonymous across the different frequency classes
```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

mkdir -p ${OUT_DIR}/HeteroticGroups

# plot stats separately for HA and RHA
group=HA
genotypes=$(awk -v var="$group" -F',' '{if ($11==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO | paste -sd,)
bcftools view -Ou --samples ${genotypes} ${VCF} |\
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${OUT_DIR}/HeteroticGroups/HA_SNP_info.txt

awk '{if ($5/$6 > 0.9) {print $0}}' ${OUT_DIR}/HeteroticGroups/HA_SNP_info.txt | wc -l # 34,593

group=RHA
genotypes=$(awk -v var="$group" -F',' '{if ($11==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_INFO | paste -sd,)
bcftools view -Ou --samples ${genotypes} ${VCF} |\
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\t%AF\n' > ${OUT_DIR}/HeteroticGroups/RHA_SNP_info.txt

awk '{if ($5/$6 > 0.9) {print $0}}' ${OUT_DIR}/HeteroticGroups/RHA_SNP_info.txt | wc -l # 200,499

# reduce SNP table to just deleterious and synonymous to save memory in R
awk '{print $13}' $SNP_INFO | sort -u
awk '{if ($13=="AllDel" || $13=="Synonymous" || $13=="Variant_type") {print $0}}' $SNP_INFO > ${OUT_DIR}/HeteroticGroups/SNP_Info_Reduced.txt

# combine 
Rscript --verbose "${REPO_DIR}/BAD_Mutations/Variant_analyses/Scripts/SNP_Freq_Groups.R" \
"${OUT_DIR}/HeteroticGroups/HA_SNP_info.txt" \
"HA" \
"${OUT_DIR}/HeteroticGroups/RHA_SNP_info.txt" \
"RHA" \
"${OUT_DIR}/IntermediateFiles/SNPINFO_ForUnfolded.txt" \
"${OUT_DIR}/HeteroticGroups/HA_RHA_DerivedFreqs.txt"
```

----- 

Continued with `dSNP_DerAncBINS.R` <- for binning across genome

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
