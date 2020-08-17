# SNP-calling Steps

All SNP-calling was performed using sequence_handling: 

### Navigation: Jump to Section

- 

--

## Re-sequencing data

Information about raw sequence data & file wrangling in: FileInformation.md  

Fastq.gz files needed to be concatenated (in some cases more than 1 forward/reverse) - see DataProcessing/Concatenate.sh

288 lines split into 9 groups of samples for processing (#1-7, "S_African_seqs", "SRA_seqs")

## Adapter Trimming

Adapter trimming with Scythe, prior of 0.05

Adapters used -
A subset of sequences show contamination with nextera transposae sequence (S_African_seqs), so trimming was redone using - 

After adapter trimming, quality assessment was re-run on trimmed samples to verify there was no residual adapter contamination

## Read Mapping

Used BWA v.0.7.17, default parameters

Read mapping statistics-

## SAM Processing

## Haplotype Caller

## GenomicsDB Import

## Genotype GVCFs

## Create HC Subset

## Variant Recalibrator

## Variant Filtering

## Variant Analysis