### My "Custom" Filtering
1.) filtered out sites that didn't pass variant recalibrator and non-variant sites (and selected only snps)
	- 81431704 out of 87332695 (93% kept)

2.) Second, filtered based on genotype fields:
	- GQ values less than 6 (10th percentile)
	- DP more than 50 (99th percentile is 25 but based on uneven coverage among samples as observed in sequence coverage graph, used a higher number)
	- Number of sites: 81,431,704 (has not changed as expected- sites just marked as filtered)
	- After filter flags in place, filtered for no more than 0.2 genotypes marked as filtered or no-call and set filtered genotypes to no-call
		-chromosome 1 before and after filtering: 4386393 v. 3245650
	- 63,739,304 out of 81431704 (78% kept)

2b.) Evaluate for QUAL and % heterozygous metrics
	- QUAL vs. ts/tv statistics
	- Use vcftools --het to calculate a measure of heterozygosity on a per-individual basis
	- Distribution of ExcessHet

### QC Statistics:
```bash
# run interactive job
srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l # Job 727860
module load BCFtools/1.10.2-GCC-8.3.0\
VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter2_122828/Sunflower_SAM_SNP_Calling_GenoFieldFiltered.vcf"
bcftools stats $VCF > SAM_SNPs_GTfilteredStats.txt
```

- QUAL vs. ts/tv
using file SAM_SNPs_filtered_Stats.txt from bcftools stat function
```bash
grep "^QUAL" SAM_SNPs_GTfilteredStats.txt > SAM_SNPs_GTfilteredStatsQUAL.txt
```

```R
QUAL <- read.table("SAM_SNPs_GTfilteredStatsQUAL.txt", header = FALSE, sep = "\t")
colnames(QUAL)[3:6] <- c("Quality", "Num_SNPs", "Num_ts", "Num_tv")
QUAL$tstv <- QUAL$Num_ts / QUAL$Num_tv
plot(QUAL$tstv ~ QUAL$Quality, cex=0.5)
abline(v=40, lty=2)
abline(h=1.71)

```

3.) Filter for QUAL > 40

- 61,182,383 out of 63,739,304 (96%); ts/tv=1.72

QC Stats
```bash
# run interactive job
srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l # Job 728964
module load BCFtools/1.10.2-GCC-8.3.0\
VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter3_123120/Sunflower_SAM_SNP_Calling_QUALFiltered.vcf"
bcftools stats $VCF > SAM_SNPs_QUALfilteredStats.txt
```
Inbreeding coefficient (F) for all individuals

Highly heterozygous individuals:
```bash
module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
vcftools --vcf $VCF --het --out SAM_hetIND_GTfiltered

module load R/3.5.0-foss-2019b
R
```
```R
het <- read.table("SAM_hetIND_GTfiltered.het", header = TRUE, sep = "\t")

het[which(het$F < 0.5),] # 35 genotypes

```

4.) Filter out heterozygous sites

This will be done in multiple steps.
First, I will 'filter' variants to mark the genotypes that are heterozygous at each site. ~~I will also use the flag "invalidate previous filters" (for the next step).~~ <- this doesn't appear to work. Will need to start from "Filter 1" (before applying the genotype quality/depth filters)

Then, I will 'select' the sites that have less than 20% of genotypes marked as heterozygous. I will then perform select Variants on the output from step 3 to select concordant sites (that aren't in the list of highly heterozygous sites).
- the 20% was selected based on: 2 rounds of inbreeding so no more than 12.5% of samples expected to be heterozygous at any given locus, however ~ 10% samples are more heterozygous than expected. 

Used `Het_filter.sh`

Number of variatns after filtering: (`Sunflower_SAM_SNP_Calling_HetFieldFiltered.vcf`): 
81,431,704 (did not filter any variants...?). Did this not work or were there no sites with >20% heterozygotes?

5.) Biallelic, remove highly heterozygous individuals