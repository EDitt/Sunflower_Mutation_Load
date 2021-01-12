### My "Custom" Filtering
1.) filtered out sites that didn't pass variant recalibrator and non-variant sites (and selected only snps)
	- 81431704 out of 87332695 (93% kept)

Used `gatk_SelectVariants.sh`

2.) Second, filtered based on genotype fields:
	- GQ values less than 6 (10th percentile)
	- DP more than 50 (99th percentile is 25 but based on uneven coverage among samples as observed in sequence coverage graph, used a higher number)
	- Number of sites: 81,431,704 (has not changed as expected- sites just marked as filtered)
	- After filter flags in place, filtered for no more than 0.2 genotypes marked as filtered or no-call and set filtered genotypes to no-call
		-chromosome 1 before and after filtering: 4386393 v. 3245650
	- 63,739,304 out of 81431704 (78% kept)

Used `gatk_FilterVarGeno.sh` and `gatk_SelectVarGeno.sh`

2b.) Evaluate for QUAL and % heterozygous metrics
	- QUAL vs. ts/tv statistics
	- Use vcftools --het to calculate a measure of heterozygosity on a per-individual basis
	- Distribution of ExcessHet

### QC Statistics:
```bash
# run interactive job
srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l # Job 727860
module load BCFtools/1.10.2-GCC-8.3.0
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

Used `gatk_SelectVariantsQUAL.sh`

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

4.) Filter out highly heterozygous sites

Use data obtained from GATK's VariantsToTable to count the number of heterozygous genotypes for each variant
```bash
srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l # Job 766736
module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk
INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120/Sunflower_SAM_SNP_Calling_snps.filtered.vcf
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121
gatk --java-options "-Xmx2g" VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F TYPE -F DP -F ExcessHet -F InbreedingCoeff -F HET -F HOM-REF -F HOM-VAR -F NCALLED \
     -O "${OUTPUT_DIR}/Variants_HetInfo.table"
```

Using GATK to filter out sites with a high proportion of heterozygotes did not work (`Het_filter.sh`)

Instead I used bcftools
Here, I'm testing it
```bash
tmux new -s bcftools_test
srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l

INPUT_VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121/Sunflower_SAM_SNP_Calling_HetFieldFiltered.vcf"
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121/test_bcftools"

module load BCFtools/1.10.2-GCC-8.3.0

bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.2' $INPUT_VCF -o ${OUTPUT_DIR}/bcftools_Het_filtered2.vcf

grep -v "^#" bcftools_Het_filtered2.vcf | wc -l # 78,629,492
grep -v "^#" $INPUT_VCF | wc -l # 81,431,704
```

I checked to make sure it did this properly using the Variants-to-Table function again from GATK.

```bash
module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk

INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121/test_bcftools/bcftools_Het_filtered2.vcf
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Het_Filter_010121/test_bcftools

# need to index first
gatk --java-options "-Xmx2g" IndexFeatureFile \
-F $INPUT_VCF

gatk --java-options "-Xmx2g" VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F TYPE -F DP -F ExcessHet -F InbreedingCoeff -F HET -F HOM-REF -F HOM-VAR -F NCALLED \
     -O "${OUTPUT_DIR}/bcftools_Variants_HetInfo.table"
```

Since this worked, I used script `Het_filter_bcftools.sh` on output from step 3

This filtered out 3.4% of sites. 
58,390,842 variants, ts/tv 1.72

5.) Will also filter out variants with high ExcessHet (>5). This number was chosen based on looking at distribution of variants in a table (see `Heterozygosity_exploration.R`).

These are sites with a lot of heterozygotes, but not many homozygous variant sites

Used script `gatk_SelectVariantsExcessHet.sh`

After filtering:
54,145,780 variants, ts/tv 1.72

6.) Select only biallelic, make sure non-variant sites removed

First, filter out non-variant and filtered sites (just a precautionary measure)
Used `gatk_SelectVarRemoveNonVar.sh`


7.) Remove highly heterozygous individuals?
