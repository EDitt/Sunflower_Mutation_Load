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

## Genotype GVCFs

## Variant Recalibrator

The Unfiltered DP per sample and GQ percentile tables are in "VCF_Stats/" as well as tables for the filtered sets
See "HC_Variant_QC" folder for VCF exploration info and information about obtaining the high confidence subsets

Training sets included:  
1.) High Confidence subset (see HC_Variant_QC folder for raw and filtered VCF statistics)  
2.) Truth Set of SNPs obtained from an Illumina SNP array (see "TruthSNPs" folder)

## Variant Filtering

## Variant Analysis