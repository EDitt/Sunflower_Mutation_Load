# SNP-calling Steps

All SNP-calling was performed using sequence_handling: 

### Navigation: Jump to Section

- [Pre-processing](#pre-processing)
- [Adapter Trimming](#adapter-trimming)
- [Read Mapping](#read-mapping)
- [SAM Processing](#sam-processing)
- [Haplotype Caller](#haplotype-caller)
- [Genomics DB Import](#genomicsdb-import)
- [Genotype GVCFs](#genotype-gvcfs)
- [Create HC Subset](#create-hc-subset)
- [Variant Recalibrator](#variant-recalibrator)
- [Variant Filtering](#variant-filtering)
---

## Pre-processing

Information about raw sequence data & file wrangling in: DataProcessing/FileInformation.md  

Many samples had more than 1 forward/reverse fastq.gz file that needed to be concatenated before beginning sequence handling - see DataProcessing/Concatenate.sh

288 lines split into 9 groups of samples for processing (#1-7, "S_African_seqs", "SRA_seqs")

---

## Adapter Trimming

Adapter trimming with Scythe, prior of 0.05

Adapters used -
A subset of sequences show contamination with nextera transposae sequence (S_African_seqs), so trimming was redone using - 

Two samples had different quality encoding (Illumina 1.5 instead of Sanger/Illumina 1.9). The quality encoding was accounted for in adapter trimming and subsequent steps

After adapter trimming, quality assessment was re-run on trimmed samples to verify there was no residual adapter contamination

---

## Read Mapping

Genome downloaded from: https://sunflowergenome.org/assembly-data/assets/data/assemblies/Ha412HOv2.0-20181130.fasta.gz
Annotations downloaded from https://www.heliagene.org/ICSG/ (Ha412-HOv2.0-20181130)
```bash
GenomeFASTA=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta.gz

gzip -dc $GenomeFASTA > /scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta
module load SAMtools/1.10-iccifort-2019.5.281 # (previously used SAMtools/1.3.1-foss-2016b)
samtools faidx /scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta
```

Used BWA v.0.7.17, default parameters

Read mapping statistics-

---

## SAM Processing

SAM Processing was performed using Picard (v.)

---

## Haplotype Caller

---

## GenomicsDB Import

In order for Genotype GVCFs to work, this step needed to be run with genomic intervals. I wrote the "Intervals_at_Ns.sh" script (sequence_handling/HelperScripts) to break up genome at strings of N's (at least 3k bp in length) into 165 regions of approx 20k. See ResourceFiles/INTERVALS_20k_atNs.bed

I also ran on Scaffold Sequences in 2 separate groups:
1.) Larger than 100k bp (N=53)
2.) 10k-100k bp (N=520) - broke into 2 separate jobs: "Medium A" and "Medium B"
  - part A list: 00054-00253 (N=200)
  - part B list: 00254-00573 (N=320)

Did not run on scaffold sequence smaller than 10kb (N=24,500)

```bash
head -70 Full_Intervals.list | tail -53 > LargeScaffolds_Intervals.list
tail -25020 Full_Intervals.list | head -520 > Med_Scaffolds_Intervals.list
```

The chromosomal regions were parallelized, while the scaffolds were run without parallelization as the 3 separate groups

---

## Genotype GVCFs

These were run on exactly the same regions as GenomicsDB Import. Scaffold sequence was run separately from chromosomal regions.

At the end, I made lists of the VCF chromosomal parts (Run in parallel) and the 3 groups of scaffolds (run as 3 separate jobs), and appended them into one list:
(I also had to change the name of the scaffold .vcf files)

```bash
Chrom=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Genotype_GVCFs/vcf_split_regions/VCF_chrom_parts_list.txt
Scaffolds=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Scaffold_Redo/ScaffoldVCF.txt

cat $Chrom $Scaffolds > VCF_parts_list_new.txt
```

---

## Create HC Subset

(When I first ran this handler, it initially filtered out indels before it did all subsequent steps. The results from this are in the "HC_Variant_QC/SNP_Only" folder)

I then ran the handler after changes were made (Sept 2020) that did all steps *with* indels included. This file was too large to get DP_per_sample files, but I used the numbers from the data from SNPs only

There were several steps to this handler. 
1.) Concatenated the split VCF file parts from Genotype GVCFs
This handler combined the VCF parts files. 


I compressed this raw file to save in `jmblab/Projct`
(this raw file also includes indels)
```bash
qsub -I -q s_interq -l walltime=24:00:00 -l nodes=1:ppn=4 -l mem=8gb

module load BCFtools/1.10.2-GCC-8.3.0
bgzip -c Sunflower_SAM_SNP_Calling_raw_variants.vcf > SAM_AlignedHA412HOv2_raw_variants.vcf.gz
```

2.) Get GQ and DP per sample percentiles
There was not enough memory to get DP_per_sample percentiles for the raw/filtered VCFs that included indels. 
I changed the following code in the `percentiles.R` script in sequence_handling to get this to run:
- line 78: changed DP_file_size limit to 50798691840
- line 81: changed print statement, "DP file is larger than 50 gb..."

3.) Filter SNPs
##### Cut-offs:
- Depth per Sample cutoff = 5
- GQ Cutoff = 6 (10th percentile of raw GQ table)
- Maximum proportion of heterozygous per site = 0.15
- Maximum proportion of missing calls (or low GQ genotypes) per site = 0.2
- QUAL score cutoff = 40


Number of SNPs/Indels in raw VCF file compared to filtered?
```bash
grep -v "#" Sunflower_SAM_SNP_Calling_raw_variants.vcf | wc -l # 101,509,931

# filtered SNPs. (Errored out at step 5 so in Intermediates directory)
grep -v "#" "Intermediates/Sunflower_SAM_SNP_Calling_filtered.vcf" | wc -l #71,298,676
 # note: there was a problem with the way the handler filtered. Will re-run after changes made

# how many SNPs versus indels?

```

Error with GATK variant filtration:
WARN  JexlEngine - ![0,2]: 'GQ < 6;' undefined variable GQ
```bash
grep "undefined variable GQ" Sunflower_SAM_SNP_Calling_Create_HC_Subset.e3336389 | wc -l # 101,509,931
```

Look at distribution of annotations

```bash
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8

gatk VariantsToTable \
     -V Sunflower_SAM_SNP_Calling_raw_variants.vcf \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F GQ \
     -O RawVariants.table
```

Used data in 'Variant_annotations.R' script

Second try (after fixing issue with handler- GATK wasn't actually filtering)
```bash
# intermediate file (after filtering for depth and quality with bcftools)
grep -v "#" Sunflower_SAM_SNP_Calling_filtered_dp_and_qual.vcf | wc -l # 97,465,016

```

Working, but not getting past step 6 due to errors with VCFtools:
`vcftools: /lib64/libstdc++.so.6: version 'GLIBCXX_3.4.21' not found (required by vcftools)`

Ran as it's own job- "HC_Subset_step6.sh"
`After filtering, kept 58054334 out of a possible 58054411 Sites`
```bash
# intermediate file (after filtering for depth and quality with bcftools)
grep -v "#" Sunflower_SAM_SNP_Calling_high_confidence_subset.vcf | wc -l # 58,054,334
```

---

## Variant Recalibrator

The Unfiltered DP per sample and GQ percentile tables are in "VCF_Stats/" as well as tables for the filtered sets
See "HC_Variant_QC" folder for VCF exploration info and information about obtaining the high confidence subsets

Training sets included:  
1.) High Confidence subset (see HC_Variant_QC folder for raw and filtered VCF statistics)  
2.) Truth Set of SNPs obtained from an Illumina SNP array (see "TruthSNPs" folder)

Output messages-

After Variant recalibrator was begun:
```bash
ProgressMeter - Traversal complete. Processed 87332695 total variants in 384.0 minutes.
INFO  VariantDataManager - QD:      mean = 19.82    standard deviation = 11.58
INFO  VariantDataManager - FS:      mean = 8.24     standard deviation = 18.99
INFO  VariantDataManager - ReadPosRankSum:          mean = 0.14     standard deviation = 0.86
INFO  VariantDataManager - MQ:      mean = 48.68    standard deviation = 8.97
INFO  VariantDataManager - MQRankSum:       mean = -0.93    standard deviation = 1.18
INFO  VariantDataManager - SOR:     mean = 1.60     standard deviation = 1.43
INFO  VariantDataManager - DP:      mean = 2379.86  standard deviation = 1327.71
INFO  VariantDataManager - Annotation order is: [DP, MQ, QD, FS, SOR, MQRankSum, ReadPosRankSum]
INFO  VariantDataManager - Training with 49808820 variants after standard deviation thresholding.
WARN  VariantDataManager - WARNING: Very large training set detected. Downsampling to 2500000 training variants.
INFO  GaussianMixtureModel - Initializing model with 100 k-means iterations...
...
INFO  VariantRecalibratorEngine - Convergence after 75 iterations!
INFO  VariantRecalibratorEngine - Evaluating full set of 87332695 variants...
INFO  VariantDataManager - Selected worst 11871283 scoring variants --> variants with LOD <= -5.0000
INFO  VariantRecalibratorEngine - Finished iteration 0.
INFO  VariantRecalibratorEngine - Finished iteration 5.    Current change in mixture coefficients = 0.14029
INFO  VariantRecalibratorEngine - Finished iteration 10.   Current change in mixture coefficients = 0.03250
INFO  VariantRecalibratorEngine - Convergence after 13 iterations!
INFO  VariantRecalibratorEngine - Evaluating full set of 87332695 variants...
INFO  TrancheManager - Finding 12 tranches for 87332695 variants

```

Exit status 127 even though it finished all steps-
At end:
`/var/spool/torque/mom_priv/jobs/3399078.sapelo2.SC: line 2: 0: command not found`
(this was confirmed by GACRC to be a bug that wasn't affecting the job status)

## Variant Filtering

This part of the handler has not been updated yet, so I will modify. 

Using GATK's SelectVariants, (`gatk_SelectVariants.sh`), I will:
1. Filter out indels
2. Excluded sites marked as filtered by ApplyVQSR
3. Exclude non-variant sites

After getting rid of indels and filtered sites, 81,431,704 variants remain (file in directory `Filter1_102120`)


#### Testing out how much remaining filtering I want to do-
1.) minimal approach (above)
	Number of records: 81,431,704
	ts/tv: 1.71
2.) Using sequence handling's Variant Filtration. 
	Default settings except
	- 1.0 for max % deviation allowed in heterozygotes (not filtering based on this parameter due to sample pooling)
	- Max. proportion of heterozygote genotypes = 0.2 (accounting for ~ 10% of lines that highly heterozygous)
3.) Some filtering above minimum- max # of reads, max number of samples un-called or low GQ, 
	minimum QUAL: 40

### Sequence handling
- After step 1: 81,431,704
- After step 3: 76,700,323
- After step 4: 65,208,272
- After step 6: 116
- After step 7: 116

### My code:



