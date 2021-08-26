# Number of dSNPs per genotype

### Setup
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

```

# Count Number of Alt/Ref Alleles per genotype
(using bcftools)
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_Refdeleterious.vcf > RefDeleterious_StatsperSample.txt
grep "PSC" RefDeleterious_StatsperSample.txt > RefDeleterious_perSampleCounts.txt

bcftools stats -s - SAM_Altdeleterious.vcf > AltDeleterious_StatsperSample.txt
grep "PSC" AltDeleterious_StatsperSample.txt > AltDeleterious_perSampleCounts.txt
```

# Add Deleterious Mutations from both lists for each genotype
```bash
module load R/4.0.0-foss-2019b
R
```
```R
DelRef_data <- read.table("RefDeleterious_perSampleCounts.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(DelRef_data) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")

#DelRef_data <- subset(DelRef_data, select = -c(PSC, id, nIndels, nHapRef, nHapAlt))

DelRef_data$Total <- DelRef_data$nRefHom + DelRef_data$nNonRefHom + DelRef_data$nHets + DelRef_data$nMissing

DelAlt_data <- read.table("AltDeleterious_perSampleCounts.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
colnames(DelAlt_data) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")

#DelAlt_data <- subset(DelAlt_data, select = -c(PSC, id, nIndels, nHapRef, nHapAlt))

DelAlt_data$Total <- DelAlt_data$nRefHom + DelAlt_data$nNonRefHom + DelAlt_data$nHets + DelAlt_data$nMissing

# merge datasets
Deldata <- merge(DelRef_data[,c(3:6,14:15)], 
					DelAlt_data[,c(3:6,14:15)],
					by="sample",
					suffixes=c("_DelRef", "_DelAlt"))
length(Deldata$sample) #288

# Number of Deleterious Alleles:
Deldata$TotNum_dSNPs <- 2*Deldata$nRefHom_DelRef +
						Deldata$nHets_DelRef +
						2*Deldata$nNonRefHom_DelAlt +
						Deldata$nHets_DelAlt

# Number of homozygous deleterious sites:
Deldata$nHom_dSNPs <- Deldata$nRefHom_DelRef +
					Deldata$nNonRefHom_DelAlt

# Number of heterozygous deleterious sites:
Deldata$nHet_dSNPs <- Deldata$nHets_DelRef +
					Deldata$nHets_DelAlt

# total number of SNPs (for proportions?)
Deldata$Total <- Deldata$Total_DelRef + Deldata$Total_DelAlt
Deldata$TotalMissing <- Deldata$nMissing_DelRef + Deldata$nMissing_DelAlt

# save file
write.table(Deldata[,c(1,12:16)], "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Genotype_dSNP_counts.txt", sep = "\t", quote=FALSE, row.names=FALSE)
```
