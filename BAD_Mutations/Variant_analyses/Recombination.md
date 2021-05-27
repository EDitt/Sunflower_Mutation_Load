# Recombination across genome

Using John Bowers genetic map information, I used the file I manipulated for the truth SNPs: `SNP_Genetic_Map_Unique.txt` for the cM distance of the markers. This file contains 6984 markers that mapped uniquely. 
- The first column is the linkage group (chromosome number), second column is the locus name, third column is the distance in cM.

I also used the vcf file that I created with SNPutils- `MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf` where I mapped the SNPs to the new genome build. (N=6523)
- Here the first column is the chromosome, 2nd is the position in bp, third is Locus name


## Create Recombination Data Files
```bash
GeneticMap=/scratch/eld72413/SNParray/SNP_Genetic_Map_Unique.txt
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

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load R/4.0.0-foss-2019b
R
```

Use R to merge datasets
```R
cM_pos <- read.table("SNParray_cMpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(cM_pos) <- c("Chromosome", "Locus", "cM")
length(cM_pos$Locus) # 6984

bp_pos <- read.table("SNParray_BPpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(bp_pos) <- c("Chromosome", "bp", "Locus")
length(bp_pos$Locus) # 6524

length(which(cM_pos$Chromosome %in% bp_pos$Chromosome))
length(which(cM_pos$Locus %in% bp_pos$Locus)) #6524
## there are spaces in the Locus names
cM_pos$Locus <- gsub(" ","",cM_pos$Locus)
bp_pos$Locus <- gsub(" ","",bp_pos$Locus)

## also spaces in the Chromosome names
cM_pos$Chromosome <- gsub(" ","",cM_pos$Chromosome)
bp_pos$Chromosome <- gsub(" ","",bp_pos$Chromosome)

Recombination <- merge(cM_pos, bp_pos, by=c("Chromosome", "Locus")) # some markers are on the scaffolds (not in chromosomes 1-17) so not included
length(Recombination$Locus) #5958

Recombination$Mbp <- Recombination$bp/1000000
Recombination$cM_Mbp <- Recombination$cM / Recombination$Mbp

# remove markers with very high cM/Mbp (>20) indicating a misalignment (previously 5958 markers, now 5950)
Recombination <- Recombination[-which(Recombination$cM_Mbp > 20),]
write.table(Recombination, file="Recombination.txt", row.names=FALSE, quote=FALSE, sep="\t")
#save(Recombination, file="Recombination.RData")
```

### Bins for average recombination across genome
```bash
### Recombination bins
cd /scratch/eld72413/SAM_seq/Recombination
# Convert Recombination.txt to a format bedtools can recognize
awk 'NR>1 {print $1"\t"$4-1"\t"$4"\t"$6}' Recombination.txt > Recombination.bed

### Number, Average, and standard deviation of recombination rate per window
bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -wo -a - -b Recombination.bed \
| bedtools groupby -i stdin \
-g 1,2,3 -c 7 -o count > Number_Markers.txt

bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -wo -a - -b Recombination.bed \
| bedtools groupby -i stdin \
-g 1,2,3 -c 7 -o mean > Ave_Recombination.txt

bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -wo -a - -b Recombination.bed \
| bedtools groupby -i stdin \
-g 1,2,3 -c 7 -o stdev > Stdev_Recombination.txt
```

Use R to merge Recombination info
I redid this after removing the 8 markers with high cM/Mbp values (>20, N=8). It only changed the mean for 1 interval- Chromosome 12 0-10Mbp (which also other markers with fairly high values)
```R
RecombBin_Num <- read.table("Number_Markers.txt",  sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
RecombBin_Ave <- read.table("Ave_Recombination.txt",  sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
RecombBin_Stdev <- read.table("Stdev_Recombination.txt",  sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)

RecombinationInfo <- merge(RecombBin_Num, merge(RecombBin_Ave, RecombBin_Stdev, by=c("V1", "V2", "V3")), 
  by=c("V1", "V2", "V3"))
colnames(RecombinationInfo) <- c("Chromosome", "Interval_start", "Interval_end", "Number_Markers",
  "Mean_cM_Mbp", "Stdev_cM_Mbp")
write.table(RecombinationInfo, file="RecombinationBins.txt", row.names=FALSE, quote=FALSE, sep="\t") # length=302
```

## Binning SNP classes across Genome

Make windows to bin count of each class of SNPs
```bash
module load BEDTools/2.29.2-GCC-8.3.0

# make genome file
awk -v OFS='\t' {'print $1,$2'} "/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta.fai" | head -17 > "/scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt"

cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
mkdir GenomicBins

### dSNPs per window
bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b SAM_deleterious.vcf -c > ./GenomicBins/Deleterious_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b SAM_tolerated.vcf -c > ./GenomicBins/Tolerated_10MbCounts.txt
#| bedtools intersect -a - -b SAM_tolerated_polarized.vcf -c > DerivedTolerated_10MbCounts.txt

bedtools makewindows -g /scratch/eld72413/SAM_seq/Recombination/GenomeFile.txt -w 10000000 \
| bedtools intersect -a - -b SAM_synonymous.vcf -c > ./GenomicBins/Synonymous_10MbCounts.txt
```

Use R to merge bin info across the 3 classes
```R
# module load R/4.0.0-foss-2019b
# R
# combine SNP info
delSNP_bins <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/Deleterious_10MbCounts.txt",  sep = "\t", header=FALSE, stringsAsFactors = FALSE)
tolSNP_bins <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/Tolerated_10MbCounts.txt",  sep = "\t", header=FALSE, stringsAsFactors = FALSE)
synonSNP_bins <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/Synonymous_10MbCounts.txt",  sep = "\t", header=FALSE, stringsAsFactors = FALSE)

All_SNP_bins <- merge(synonSNP_bins, merge(tolSNP_bins, delSNP_bins, by=c("V1", "V2", "V3")), by=c("V1", "V2", "V3"))
colnames(All_SNP_bins) <- c("Chromosome", "Interval_start", "Interval_end", "NumSynon", "NumTol", "NumDel")
length(All_SNP_bins$Chromosome) #325

# write file
write.table(All_SNP_bins, file="/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/GenomicBins/SNP_bins", row.names=FALSE, quote=FALSE, sep="\t") 

```



On local computer:

```{r code}


# long format
Recomb_dSNP_long <- gather(Recomb_dSNP, Type, Number, PropdSNP:PropTol, factor_key=TRUE)

p <- ggplot(data=Recomb_dSNP_long, aes(x=Mean_cM_Mbp, y=Number, group=Type))
p + geom_point(aes(col=Type)) + geom_smooth(method='lm', aes(col=Type))

q <- ggplot(data=Recomb_dSNP[which(Recomb_dSNP$Mean_cM_Mbp < 4),], aes(x=Mean_cM_Mbp, y=PropdSNP_total))
q + geom_point() + geom_smooth(method='lm')

mod1 <- lm(PropdSNP_total ~ Mean_cM_Mbp +
Chromosome, data=Recomb_dSNP[which(Recomb_dSNP$Mean_cM_Mbp < 4),])
summary(mod1) # -0.004868, p=0.02225
Anova(mod1, test = "Chi")

plot(PropdSNP_total ~ Mean_cM_Mbp, 
data=Recomb_dSNP[which(Recomb_dSNP$Mean_cM_Mbp < 4),])
abline(mod1)

### greater than 10 markers represented?
hist(Recomb_dSNP$Number_Markers)

mod2 <- lm(PropdSNP_total ~ Mean_cM_Mbp +
Chromosome, data=Recomb_dSNP[which(Recomb_dSNP$Mean_cM_Mbp < 4 &
Recomb_dSNP$Number_Markers > 10),])
summary(mod2) # -0.002, p<0.0001
Anova(mod2, test = "Chi") #ns

plot(PropdSNP_total ~ Mean_cM_Mbp, 
data=Recomb_dSNP[which(Recomb_dSNP$Mean_cM_Mbp < 4 &
Recomb_dSNP$Number_Markers > 10),])
abline(mod2)
```

Plot across genome
```{r code}


plot(Recombination_Chroms$Ha412HOChr01$Position, 
      Recombination_Chroms$Ha412HOChr01$Mean_cM_Mbp, type="l")

p <- ggplot(data=Recombination_Chroms$Ha412HOChr01, aes(x=Position, y=Mean_cM_Mbp))
p1 <- p + geom_line(col="black") + 
          geom_line(aes(y=Mean_cM_Mbp - SE), col="darkgrey", linetype = "dashed") +
          geom_line(aes(y=Mean_cM_Mbp + SE), col="darkgrey", linetype = "dashed") +
          theme_minimal() + ylab("Average cM/Mbp") + xlab("")
    #geom_errorbar(aes(ymin=Mean_cM_Mbp - SE, ymax=Mean_cM_Mbp + SE),      
    #position=position_dodge(200000000), width=2000000) +

p2 <- ggplot(data=Recombination_Chroms$Ha412HOChr01, aes(x=Position, y=PropdSNP_total))
p1b <- p2 + geom_line(col="red") +
               theme_minimal() + ylab("Proportion dSNPs") + xlab("Position (Mb)")
               #+ ylim(0,0.12)

test1 <- ggarrange(p1, p1b, ncol=1, nrow=2)

```


## Older Code

Use R to merge bin info across the 3 classes w/ recombination info <- this file was done with the polarized SNP bins, later changed
```R
# combine with dSNP info
Derived_dSNP_bins <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/DerivedDeleterious_10MbCounts.txt",  sep = "\t", header=FALSE, stringsAsFactors = FALSE)
Derived_tolSNP_bins <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/DerivedTolerated_10MbCounts.txt",  sep = "\t", header=FALSE, stringsAsFactors = FALSE)

both_SNP_bins <- merge(Derived_tolSNP_bins, Derived_dSNP_bins, by=c("V1", "V2", "V3"))
colnames(both_SNP_bins) <- c("Chromosome", "Interval_start", "Interval_end", "NumTol", "NumDel")
length(both_SNP_bins$Chromosome) #325

# merge with recombination
Recombination_dSNP_bins <- merge(RecombinationInfo, both_SNP_bins, by=c("Chromosome", "Interval_start", "Interval_end"),
  all=TRUE)
length(Recombination_dSNP_bins$Chromosome) #325

# write file
write.table(Recombination_dSNP_bins, file="/scratch/eld72413/SAM_seq/dSNP_results/Recombination_dSNP_bins", row.names=FALSE, quote=FALSE, sep="\t") 

```