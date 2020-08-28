## Using ANGSD-wrapper to run Pop Gen analyses: https://github.com/ANGSD-wrapper/angsd-wrapper  

#### Pre-processing

Sequence data was aligned to the HA412v.2 genome - see SNP-callling for H. annuus data and Outgroups for other Helianthus spp.  
All data were processed using the "SAM_Processing" handler from sequence handling.  
The resulting .bam files were then used in 'Realigner_Target_Creator' and 'Indel_Realigner' using sequence handling and GATK v. 3.8-1

#### Setup

To install dependencies, first needed to load gsl module
```bash
module load GSL/2.6-GCC-8.3.0

# other dependencies needed:
module load SAMtools/1.10-GCC-8.3.0
module load gnuplot/5.2.2-foss-2018a
```

Interval List
Format interval list for chromosomal sequence for ANGSD]
```bash
INTERVALS20k=/scratch/eld72413/SAM_seq/results2/VCF_results_new/N_Intervals/INTERVALS_20k_atNs.bed

awk '{print $1":"$2"-"$3}' $INTERVALS20k > Chromsome_regions.txt

# for testing (first 5 chromosomes):
head -49 Chromsome_regions.txt > Chromosome_regionsTest.txt

```

Then I was able to clone the repository (latest commit `6d10630`) and install dependencies

#### 1. Ancestral Sequence

```bash
./angsd-wrapper Ancestral /home/eld72413/DelMut/Sunflower_Mutation_Load/ANGSD/ConfigFiles/Ancestral_Sequence_Config 
```

#### 2. Inbreeding
Tested with subset of samples (Group 7, N=9) using the first 5 chromosomes as a regions file
Using interactive job
```bash
qsub -I -q s_interq -l walltime=12:00:00 -l nodes=1:ppn=8 -l mem=22gb
```

```bash
./angsd-wrapper Inbreeding /scratch/eld72413/NSFproj/ANGSD_FILES/Inbreeding_Coefficients_Config
```