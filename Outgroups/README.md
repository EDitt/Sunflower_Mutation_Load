# Sunflower Outgroups

Reads were downloaded from the Sequence Read Archive from BioProject PRJNA397453

Sequence data was obtained as described in Hubner et al. (2019)

---

### Mapping sunflower outgroups  

_Helianthus debilis_ :

| SRA Name  | Renamed as | Mapped 4% | Mapped 3% | Mapped 5% |
|-----------| ---------- | ----------| ----------| ----------|
|SRS2413722 | Debilis_22 | 97.53%	 | 97.54%*	 | 97.51%	 |
|SRS2413743 | Debilis_43 | 97.83%	 | 98.05%*	 | 98.02%	 |
|SRS2413744 | Debilis_44 | 97.75%*	 | 97.69%	 | 97.65%	 |
|SRS2413741 | Debilis_41 | 97.76%	 | 97.86%*	 | 97.82%	 |
|SRS2413740 | Debilis_40 | 97.84%*	 | 97.74%	 | 97.73%	 |
|SRS2413739 | Debilis_39 | 97.67%	 | 97.85%*	 | 97.81%	 |
|SRS2413737 | Debilis_37 | 98.03%*	 | 97.77%	 | 97.73%	 |

				Average:	97.77%		97.79%		97.75%
---
* one of the 5% .sam files is truncated

## Methods

Reads were downloaded using `fasterq-dump` from the SRA toolkit (2.9.6)
> See `SRA_download.sh`

Samples were quality assessed and trimmed using sequence handling https://github.com/EDitt/sequence_handling

Samples were then mapped to the HA412Hov.2.0 genome
- First, a genome index file was created using the following commands:
```bash
module load Stampy/1.0.31-foss-2016b-Python-2.7.14 
cd /scratch/eld72413/NSFproj/ancestralseqs/GenomeFiles/
stampy.py -G Ha412HOv2 /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta  
# (The `-G` flag specifies the PREFIX for the genome index)  
# This created a genome file "Ha412HOv2.stidx"  
# Then a genome hash file was created in the same directory:  
stampy.py -g Ha412HOv2 -H Ha412HOv2
# Inputs needed are `-g` to specify the genome index prefix (use the genome index file PREFIX.stidx), and `-H` to build a hash file with the prefix listed (build hash PREFIX.sthash)
```

- Finally samples were mapped using the `Stampy_mapping.sh` script. Samples were mapped with a 0.03, 0.04, and 0.05 substitution rate to check mapping quality
- Used samtools flagstat to check mapping percentages:
```bash
module load SAMtools/1.10-GCC-8.3.0
for file in *.sam; do
	samtools flagstat $file >> Mapped.05_Debilis.txt
done
```
Mapping percentage was highest with a 0.03 substitution rate so these files were used for SAM_Processing  

SAM_Processing with Picard was done using sequence handling. Final stats:

|Accession	|	Reads_Mapped|
 ----------- ----------------
|22	|	14,137,927 |
|37	|	36,615,622 |
|39	|	34,510,330 |
|40	|	32,270,543 |
|41	|	54,225,956 |
|43	|	36,648,908 |
|44	|	38,410,001 |

Chose accession with highest number of mapped reads (#41) to use for ancestral sequence

--

## Helianthus annuus Landraces
Resequencing data for 20 H. annuus landraces is available on SRA
- Landraces already included in SAM lines: SAM046 (Mandan #1); Hopi dye (SAM083), HOPI (SAM285). Deleted from the list of 20 ("Hopi", "SAM083", "SAM046")

```bash
awk '{print $3}' ./Sunflower_Mutation_Load/Outgroups/Annuus_Landraces > Landrace_SRA.list
```

Used SRA_download.sh to obtain the sequence data for 17 landraces (not included in SAM lines)

--

## Wild Helianthus annuus

List of individuals collected in Todesco et al. 2020
4 entire populations
- Chose populations to sample from with largest number of pops and obtained SRA numbers (see SRA_Metadata.R)

- Representative Individuals - chose individual with 
```bash
awk '{print $4}' WildAnnuusPops
```


