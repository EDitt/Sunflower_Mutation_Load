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
--

## Pre-processing

Information about raw sequence data & file wrangling in: DataProcessing/FileInformation.md  

Many samples had more than 1 forward/reverse fastq.gz file that needed to be concatenated before beginning sequence handling - see DataProcessing/Concatenate.sh

288 lines split into 9 groups of samples for processing (#1-7, "S_African_seqs", "SRA_seqs")

## Adapter Trimming

Adapter trimming with Scythe, prior of 0.05

Adapters used -
A subset of sequences show contamination with nextera transposae sequence (S_African_seqs), so trimming was redone using - 

Two samples had different quality encoding (Illumina 1.5 instead of Sanger/Illumina 1.9). The quality encoding was accounted for in adapter trimming and subsequent steps

After adapter trimming, quality assessment was re-run on trimmed samples to verify there was no residual adapter contamination

## Read Mapping

Used BWA v.0.7.17, default parameters

Read mapping statistics-

## SAM Processing

SAM Processing was performed using Picard (v.)

## Haplotype Caller


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

## Genotype GVCFs


## Create HC Subset

There were several steps to this handler. First, it concatenated the split VCF file parts from Genotype GVCFs
This handler combined the VCF parts files. I compressed this raw file to save in `jmblab/Projct`
(this raw file also includes indels)

```bash
module load BCFtools/1.10.2-GCC-8.3.0
bgzip -c Sunflower_SAM_SNP_Calling_raw_variants.vcf > SAM_AlignedHA412HOv2_raw_variants.vcf.gz
```

## Variant Recalibrator

The Unfiltered DP per sample and GQ percentile tables are in "VCF_Stats/" as well as tables for the filtered sets
See "HC_Variant_QC" folder for VCF exploration info and information about obtaining the high confidence subsets

Training sets included:  
1.) High Confidence subset (see HC_Variant_QC folder for raw and filtered VCF statistics)  
2.) Truth Set of SNPs obtained from an Illumina SNP array (see "TruthSNPs" folder)

## Variant Filtering

## Variant Analysis