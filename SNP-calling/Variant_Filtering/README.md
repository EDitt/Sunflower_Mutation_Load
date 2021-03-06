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

6.) Make sure non-variant sites removed

First, filter out non-variant and filtered sites (just a precautionary measure)
Used `gatk_SelectVarRemoveNonVar.sh`

54,145,780 variants remain (did not remove any, which was expected)

* new: decided to filter for min DP per sample
Make a script `gatk_FilterMinDP.sh` to run
```bash
ERROR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921"
sbatch -o ${ERROR}/Filter.%j.out -e ${ERROR}/Filter.%j.err gatk_FilterMinDP.sh #1542394
```
45,973,312 variants remain



7.) Select only biallelic

Filter to biallelic sites
```bash
tmux new -s biallelic
srun --pty  -p inter_p  --mem=8G --nodes=1 --ntasks-per-node=4 --time=24:00:00 --job-name=qlogin /bin/bash -l #1027883

module load BCFtools/1.10.2-GCC-8.3.0
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Biallelic"
VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf"

bcftools view -m2 -M2 -v snps --threads 4 ${VCF} --output-type v --output-file ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_BIALLELIC.vcf

# this completed-  51,014,412, ts/tv ratio is 1.82
```


8.) Remove highly heterozygous individuals?


----

## All steps at once

Subset VCF to test first
```bash
module load BCFtools/1.10.2-GCC-8.3.0
out_dir=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test
# 25141 header lines
bcftools view /scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz | head -60000 > ${out_dir}/Test.vcf
# 34,859 sites
```

test on small file
```bash
ERROR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test/ErrorFiles
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test/Test.vcf' -o ${ERROR}/Filter.%j.out -e ${ERROR}/Filter.%j.err Filter_AllSteps.sh #1542038
# [filter.c:2491 filters_init1] Error: the tag "${het_prop}" is not defined in the VCF header
```

Run on recalibrated file
Will direct standard output to be saved here:
`/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/ErrorFiles`

```bash
ERROR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/ErrorFiles
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz' -o ${ERROR}/Filter.%j.out -e ${ERROR}/Filter.%j.err Filter_AllSteps.sh #1542380
# job ran out of walltime (48 hours)
# restarted after commenting out steps that had been completed- new: 1566818

# with multi-allelic sites, had 52707200 sites; biallelic= 49596706

```

Output:
Selecting Pass Sites from Variant Recalibrator
After filtering out sites that failed variant recalibrator, there are 81431704 sites left

Adding Filter labels to low quality genotypes. Cut-offs are as follows: MinGQ is 6; MinDP is 3; MaxDP is 50
Done adding filter labels to low quality genotypes. Removing sites with more than 0.2 low quality or missing genotypes
After filtering out sites with too many low quality or missing variants, there are 62245682 sites left

Filtering out sites with more than 0.2 heterozygous genotypes
After filtering out sites with more than 0.2 heterozygous sites, there are 59273737 sites remaining

Removing sites with QUAL values below 40.0 and/or ExcessHet values above 5.0
After filtering for QUAL and ExcessHet annotations, there are 52707200 sites remaining

Selecting biallelic sites
After selecting only biallelic sites, there are 49596706 variants

Why is there this discrepancy? Filtering with a different order (see step 7) gave me 45,973,312 multi-allelic sites
(compared to 52707200)

Look at differences with variants to table
```bash
#tmux window: count
module load BCFtools/1.10.2-GCC-8.3.0
bcftools stats Sunflower_SAM_SNP_FINALFilter.vcf > SAM_SNPs_FINAL_multiallelic.txt

srun --pty  -p inter_p  --mem=8G --nodes=1 --ntasks-per-node=4 --time=24:00:00 --job-name=qlogin /bin/bash -l

module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
GATK_JAR=/apps/eb/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8/gatk
INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_FINALFilter.vcf
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All
gatk --java-options "-Xmx2g" VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F FILTER -F DP -F ExcessHet -F QUAL -F InbreedingCoeff -F HET -F HOM-REF -F HOM-VAR -F NCALLED -F NO-CALL -F MULTI-ALLELIC \
     -GF GQ -GF DP \
     --show-filtered \
     -O "${OUTPUT_DIR}/Variants_VarFilterAllFINAL.table"
```

and on new filter:
```bash
INPUT_VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921/Sunflower_SAM_SNP_Calling_DP_min3Filtered.vcf
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed
OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921
gatk --java-options "-Xmx2g" VariantsToTable \
     -V "${INPUT_VCF}" \
     -L "${INTERVALS}" \
     -F CHROM -F POS -F FILTER -F DP -F ExcessHet -F QUAL -F InbreedingCoeff -F HET -F HOM-REF -F HOM-VAR -F NCALLED -F NO-CALL -F MULTI-ALLELIC \
     -GF GQ -GF DP \
     --show-filtered \
     -O "${OUTPUT_DIR}/Variants_DPmin3Filtered.table"
```

One problem I uncovered for both sets was the degree of missingness. They both have sites with a much higher number of missing data than 0.2. This might be to GATK performing that filtering step *before* it sets filtered genotypes to no-call.

I made a script - `bcftools_missingfilter.sh` to filter missing data from some of the intermediate vcf files

```bash
sbatch --export=vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Intermediates/Sunflower_SAM_SNP_GenoFieldFiltered.vcf',out_dir='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Intermediates',out_prefix='Sunflower_SAM_SNP' bcftools_missingfilter.sh # 1892317
```
Number = 44,689,089 (previously was 62,245,682)

I will also do this on the other vcf file (different order):
```bash
sbatch --export=vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921/Sunflower_SAM_SNP_Calling_DP_min3Filtered.vcf',out_dir='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter7_020921',out_prefix='Sunflower_SAM_SNP' bcftools_missingfilter.sh # 1892342
```
Number = 39,245,434 (previously was 45,973,312)

# Re-ran after adding new step 3 (bcftools) to filter out sites with too many missing variants (set to no-call by GATK):
```bash
# test first
ERROR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test/ErrorFiles
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Test/Test.vcf' -o ${ERROR}/Filter.%j.out -e ${ERROR}/Filter.%j.err Filter_AllSteps.sh # 1897190

# deleted the "Het filtered", "FINAL filter", and "Biallelic filter" vcf files and restarted from the new intermediate "missing filtered" - step 4 (redo Het filter, QUAL filter + ExcessHet, Biallelic filters)
ERROR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/ErrorFiles
sbatch --export=INPUT_VCF='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Variant_Recalibrator/Sunflower_SAM_SNP_Calling_snps.recalibrated.vcf.gz' -o ${ERROR}/Filter.%j.out -e ${ERROR}/Filter.%j.err Filter_AllSteps.sh # 1897249
```
#### Output:
Recalibrated pass sites already selected, proceeding to step 2
Filter labels already added to low quality genotypes
Filtered genotypes already set to no call, proceeding to step 3
Sites with too many low quality genotypes (set to missing by GATK in previous step) already filtered out, proceeding to step 4
Filtering out sites with more than 0.2 heterozygous genotypes
After filtering out sites with more than 0.2 heterozygous sites, there are 43034250 sites remaining
Removing sites with QUAL values below 40.0 and/or ExcessHet values above 5.0
Tool returned:
/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Intermediates/Sunflower_SAM_SNP_HETFiltered.vcf.idx
After filtering for QUAL and ExcessHet annotations, there are 39226506 sites remaining
Selecting biallelic sites
After selecting only biallelic sites, there are 37129915 variants

## Check filtering
```bash
tmux new -s VarTable
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=4:00:00 --job-name=qlogin /bin/bash -l
module load GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8

OUTPUT_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All
INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Create_HC_Subset/Intermediates/Genome_Random_Intervals.bed

gatk --java-options "-Xmx2g" IndexFeatureFile \
     -F /scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_BIALLELIC.vcf

gatk VariantsToTable \
     -V /scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_BIALLELIC.vcf \
     -F CHROM -F POS -F MULTI-ALLELIC -F QUAL -F FILTER -F NCALLED -F HET -F ExcessHet -F DP -F InbreedingCoeff -F QD \
     -GF GQ -GF DP \
      -L "${INTERVALS}" \
     --show-filtered \
     -O "${OUTPUT_DIR}/BIALLELIC_Variants.table"

# break up df to load into R
cut -f1-11 BIALLELIC_Variants.table > BIALLELIC_Variants_SITEFIELDs.table

module load R/4.0.0-foss-2019b
R

# GQ and DP values
# GQ
GQ_columns=$(head -n 1 BIALLELIC_Variants.table | tr '\t' '\n' | cat -n | grep -E *.GQ | awk '{print $1}' | paste -s -d, -)
cat BIALLELIC_Variants.table | cut -f${GQ_columns} > BIALLELIC_Variants_GQvalues.table

#DP
DP_columns=$(head -n 1 BIALLELIC_Variants.table | tr '\t' '\n' | cat -n | grep -E *.DP | awk '{print $1}' | paste -s -d, -)
cat BIALLELIC_Variants.table | cut -f${DP_columns} > BIALLELIC_Variants_DPvalues.table
```
Use R on table data
```R
NewFilter <- read.table("BIALLELIC_Variants_SITEFIELDs.table", header = T, na.strings=c("","NA"), sep = "\t")
# check filters:
# make sure no multi-allelic or failed sites
aggregate(NewFilter$POS, by=list(NewFilter$MULTI.ALLELIC, NewFilter$FILTER), length) # all "false" for multi-allelic and "PASS" for filter

# number of called genotypes
min(NewFilter$NCALLED) # 231 

# proportion of heterozygotes
NewFilter$PropHET <- NewFilter$HET / NewFilter$NCALLED
max(NewFilter$PropHET) # 0.1992883

# QUAL & ExcessHet values
min(NewFilter$QUAL) # 40.01
max(NewFilter$ExcessHet) # 4.9994

GQvals <- read.table("BIALLELIC_Variants_GQvalues.table", header = T, na.strings=c("","NA"), sep = "\t") # 2.61 GB memory
# there are values below 6 here. are those genotypes missing but GQ value still there?

GQvals[GQvals < 6] <- NA # replace with NA (assuming annotation is there even if genotype is not?)

# length of GQ values lower than the threshold
GQvals$lengthGQ <- apply(GQvals, 1, function(x) sum(is.na(x)))
# there shouldn't be more than 57
max(GQvals$lengthGQ) # 57

GQvals$MeanGQ <- rowMeans(GQvals, na.rm=T)

GQ_Info <- GQvals[,c("lengthGQ", "MeanGQ")]
colnames(GQ_Info) <- c("Num_LowQQ", "MeanGQ")

rm(GQvals) # to free up memory

DPvals <- read.table("BIALLELIC_Variants_DPvalues.table", header = T, na.strings=c("","NA"), sep = "\t")
DPvals[DPvals < 3] <- NA
DPvals$Num_LowDP <- apply(DPvals, 1, function(x) sum(is.na(x)))

DPvals[DPvals > 50] <- NA
DPvals$Num_filteredDP <- apply(DPvals, 1, function(x) sum(is.na(x))) # all filtered
max(DPvals$Num_filteredDP) # 59 ???????

DPvals$Num_HighDP <- DPvals$Num_filteredDP - DPvals$Num_LowDP

DPvals$MeanDP <- rowMeans(DPvals, na.rm=T)

### make a new dataframe

GT_summaries <- cbind(DPvals[,c("Num_LowDP", "Num_HighDP", "Num_filteredDP", "MeanDP")], GQ_Info[,c("Num_LowQQ", "MeanGQ")])

write.csv(GT_summaries, "GT_summarySTATS.csv")

### average for each accession
Ave_DP <- apply(DPvals[,c(2:289)], 2, function(x) {mean(x, na.rm=TRUE)})
write.table(Ave_DP, "AverageDP_perLine.txt")

```
Other than having 59 NA's (max should be 57) for the depth genotype fields, everything looks as expected-

## Gzip and Index VCF for VeP steps
```bash
sbatch --export=file='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_BIALLELIC.vcf' gzip_vcf.sh # 1902364
```