## Using ANGSD-wrapper to run Pop Gen analyses: https://github.com/ANGSD-wrapper/angsd-wrapper  

#### Pre-processing

Sequence data was aligned to the HA412v.2 genome - see SNP-callling for H. annuus data and Outgroups for other Helianthus spp.  
All data were processed using the "SAM_Processing" handler from sequence handling.  
The resulting .bam files were then used in 'Realigner_Target_Creator' and 'Indel_Realigner' using sequence handling and GATK v. 3.8-1

#### Setup

To install dependencies, first needed to load gsl module
```bash
module load GSL/2.6-GCC-8.3.0
module load GSL/2.6-iccifort-2019.5.281 # trying this after having complition issues

# other dependencies needed:
#module load SAMtools/1.10-GCC-8.3.0
#module load gnuplot/5.2.2-foss-2018a
module load SAMtools/1.10-iccifort-2019.5.281
module load gnuplot/5.2.8-GCCcore-8.3.0
```


##### Genomic Regions File
###### Update 01/20/21:
Will focus on gene space. Make genome file with all gene space plus 20% around the edges. 
Gene Regions + 20% of region on both sides

```bash
OUTPUTDIR="/scratch/eld72413/SAM_seq/ANGSD/Intervals"
# make a genome file:
awk -v OFS='\t' {'print $1,$2'} "/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta.fai" > "${OUTPUTDIR}/GenomeFile.txt"
# take genic regions
awk '{if ($3 == "gene") {print $0}}' /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3 > "${OUTPUTDIR}/GeneRegions.gff3"

module load BEDTools/2.29.2-GCC-8.3.0

bedtools merge -i "${OUTPUTDIR}/GeneRegions.gff3" > "${OUTPUTDIR}/GeneRegionsMerged.bed" 

bedtools slop -i "${OUTPUTDIR}/GeneRegionsMerged.bed" -g "${OUTPUTDIR}/GenomeFile.txt" -b 0.2 -pct > "${OUTPUTDIR}/GenicIntervals.bed"


# convert to regions file for angsd
awk '{print $1":"$2"-"$3}' "${OUTPUTDIR}/GenicIntervals.bed" > RegionsFile_genes.txt
```

##### Sample List
Using 288 cultivated sunflower alignment files for which indel realignment has been performed with GATK
```bash
BAMDIR="/scratch/eld72413/SAM_seq/BAM_realigned"

find ${BAMDIR} -name "*.bam"  | sort -V > "/scratch/eld72413/SAM_seq/ANGSD/SampleList.txt"

```

#### 1. Ancestral Sequence

```bash
./angsd-wrapper Ancestral /home/eld72413/DelMut/Sunflower_Mutation_Load/ANGSD/ConfigFiles/Ancestral_Sequence_Config 
```
* this was saved in my Project folder from previous analyses and I copied it into the Scratch working directory


#### 2. Inbreeding
Edited Angsd script to reflect Slurm
```bash

srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l

--chdir=/scratch/eld72413/Tmp

./angsd-wrapper Inbreeding /scratch/eld72413/SAM_seq/ANGSD/Configuration_Files/Inbreeding_Coefficients_Config
./angsd-wrapper Inbreeding /home/eld72413/DelMut/Sunflower_Mutation_Load/ANGSD/ConfigFiles/...

sbatch --export=WRAPPER='Inbreeding',CONFIG='/scratch/eld72413/SAM_seq/ANGSD/Configuration_Files/Inbreeding_Coefficients_Config' ANGSD_Job.sh
```
check processor architecture from running interactive job:
```bash
lscpu
# 2 threads per core
# tried setting N_Cores=2 in common config
```
#### 3. SFS
```bash
./angsd-wrapper SFS /scratch/eld72413/NSFproj/ANGSD_FILES/Site_Frequency_Spectrum_Config
```

### A smaller subset of genome for NSF presentation:
```bash
module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

bedtools random -l 10000 -n 150 -seed 65 -g ChromosomeLengths.txt | sort -V | awk '{print $1":"$2"-"$3}' > GenomeSubset/Random150x10k_regions.txt
```

### Combining Cultivated and Wild for Admixture analysis
```bash
CultBAMs=/scratch/eld72413/SAM_seq/BAM_realigned/Subset1_BamRealigned.txt
WildBAMs=/scratch/eld72413/NSFproj/ancestralseqs/Annuus/Indel_Realigner/Wild_RealignedBams.txt

Cult_INBREEDING=/scratch/eld72413/NSFproj/ANGSD_FILES/GenomeSubset/Cultivated_GenomeSubset/Inbreeding_Coefficients/Cultivated_GenomeSubset.indF
Wild_INBREEDING=/scratch/eld72413/NSFproj/ANGSD_FILES/GenomeSubset/Wild_GenomeSubset/Wild_GenomeSubset/Inbreeding_Coefficients/Wild_GenomeSubset.indF

cat $CultBAMs $WildBAMs > CultWildBAMs.txt
cat $Cult_INBREEDING $Wild_INBREEDING > CultWildInbreeding
```

Submitted job for Genotype likelihood analysis


## NEW ANALYSES- DEC 2020
```bash
module load GSL/2.6-iccifort-2019.5.281

# other dependencies needed:
#module load SAMtools/1.10-GCC-8.3.0
module load gnuplot/5.2.8-GCCcore-8.3.0

./angsd-wrapper setup dependencies
```


## PREVIOUS ANALYSES WITH DATA SUBSETS:

Interval List
Format interval list for chromosomal sequence for ANGSD
(below is when testing)
```bash
INTERVALS20k=/scratch/eld72413/SAM_seq/results2/VCF_results_new/N_Intervals/INTERVALS_20k_atNs.bed

awk '{print $1":"$2"-"$3}' $INTERVALS20k > Chromsome_regions.txt

# for testing (first 5 chromosomes):
head -49 Chromsome_regions.txt > Chromosome_regionsTest.txt

```

Interval list for chromosomal sequence for ANGSD
It is highly recommended to use an intervals file, as ANGSD is computationally expensive. I will use bedtools to choose random intervals for this. (later I might pipe into `bedtools shuffle` to exclude strings of N's)  
The genome file here is chromosomal sequence only
```bash
FASTA_INDEX=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta.fai

# genome file for bedtools that contains only chromosomal sequence
head -17 $FASTA_INDEX | awk -v OFS='\t' '{print $1,$2}' > ChromosomeLengths.txt

# use bedtools random
module load BEDTools/2.29.2-GCC-8.2.0-2.31.1

bedtools random -l 100000 -n 250 -seed 56 -g ChromosomeLengths.txt | sort -V | awk '{print $1":"$2"-"$3}' > Random250x100k_regions.txt
```

To look at specific chromosomal regions:
Chromosome #10 has the branching locus
```bash
awk -F "[:,-]" '{$1=$1; if ($1 == "Ha412HOChr10") {print $0}}' Random250x100k_regions.txt | awk '{print $1":"$2"-"$3}' > Chrom10/Random250x100k_Chrom10.txt

```

Then I was able to clone the repository (latest commit `6d10630`) and install dependencies

Subset 1 - 138 samples
```bash
find $(pwd -P) -name "*.bam" | sort -V > Subset1_BamRealigned.txt

# some samples needed to be re-indexed
parallel/20200422-GCCcore-8.3.0
module load SAMtools/1.10-GCC-8.3.0

find -maxdepth 1 -name "*.bam" | parallel samtools index {}


for f in `find -maxdepth 1 -name "*.bam"`; do
	echo $f
	samtools index $f
done

```
