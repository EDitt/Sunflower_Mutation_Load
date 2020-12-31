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

3.)Evaluate for QUAL and % heterozygous metrics
	- QUAL vs. ts/tv statistics
	- Use vcftools --het to calculate a measure of heterozygosity on a per-individual basis
	- Distribution of ExcessHet
- Filter ...(QUAL heterozygosity)
- Biallelic, remove highly heterozygous individuals


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
```

Highly heterozygous individuals:
```bash
module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
vcftools --vcf $VCF --het --out SAM_hetIND_GTfiltered
```