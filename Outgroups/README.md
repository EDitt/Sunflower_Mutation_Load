# Sunflower Outgroups

Reads were downloaded from the Sequence Read Archive from BioProject PRJNA397453

Sequence data was obtained as described in Hubner et al. (2019)

---

### Mapping sunflower outgroups  

_Helianthus debilis_ :

| SRA Name  | Renamed as | Mapped 4% | Mapped 3% |	
|-----------| ---------- | ----------| ----------|
|SRS2413722 | Debilis_22 | 97.53%	 |
|SRS2413743 | Debilis_43 | 97.83%	 |
|SRS2413744 | Debilis_44 | 97.75%	 |
|SRS2413741 | Debilis_41 | 97.76%	 |
|SRS2413740 | Debilis_40 | 97.84%	 |
|SRS2413739 | Debilis_39 | 97.67%	 |
|SRS2413737 | Debilis_37 | 98.03%	 |

---

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

- Finally samples were mapped using the `Stampy_mapping.sh` script. A 0.04 substitution rate was used to start with
- Used samtools flagstat to check mapping percentages:
```bash
module load SAMtools/1.10-GCC-8.3.0
for file in *.sam; do
	samtools flagstat $file >> Mapped.04_Debilis.txt
done
```

I then used SAM_Processing with sequence handling