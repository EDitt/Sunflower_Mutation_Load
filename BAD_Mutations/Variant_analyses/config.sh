#!/bin/bash

###########################################
##########     RESOURCE FILES    ##########
###########################################
#	Resource files that are common across analyses

# Final full VCF file produced by sequence handling, contains 37 M. variants
VCF=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
	# permanantly saved at: /project2/jmblab/dittmar/FinalSNPsets/Filtered

# The output of VeP run on the final, filtered VCF file
VEP=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
	# permanently saved at: /project2/jmblab/dittmar/NSF_Proj/dSNP_files

# The compiled report output from BAD_Mutations
COMPILED=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt
	# permanently saved at: /project2/jmblab/dittmar/NSF_Proj/dSNP_files

# Ancestral FASTA file created by ANGSD
ANCESTRAL=/scratch/eld72413/SAM_seq/ANGSD/Ancestral/SRS2413741_0.03_realigned.fa
	# permanantly saved at: /project2/jmblab/dittmar/ANGSD/Debilis_ancestral

# Spreadsheet with genotype lines and grouping information. "SequenceName" column is the same name that appears in the VCF
SAM_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/LineKeywINFO.csv
	# permanently saved on local computer

###########################################
########     CREATED RESOURCES    #########
###########################################
#	Files created by scripts in the current directory

# Directory with lists of SNP positions for different variant classes
POSITIONS_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/AlleleClassVCFs/FinalPositionFiles

# Ancestral alleles at all sites in VCF (only 13.5 M out of 37 M with calls). Used R to parse a .bed file
ANCESTRAL_STATE=/scratch/eld72413/SAM_seq/Polarized/AncestralStateTable
	# permanently saved at: /project2/jmblab/dittmar/NSF_Proj/SupportingData
	# only sites with calls: `awk '{if ($3 != "N") {print $0}}' ${ANCESTRAL_STATE} > AncestralStateCalls.txt`

# VEP report + Compiled preditions + Ancestral State Calls/Category
DSNP_DATA=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data_Polarized.table
	# permanently saved at: /project2/jmblab/dittmar/NSF_Proj/dSNP_files

# SNP Position, Ref/Alt alleles, Frequency Info, 
SNP_INFO=/scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/All_SNP_Info.txt
	# permanently saved at: /project2/jmblab/dittmar/NSF_Proj/dSNP_files

###########################################
#########     COMMON VARIABLES    #########
###########################################

# the location of the sunflower mutation load repository directory
REPO_DIR=/home/eld72413/DelMut/Sunflower_Mutation_Load

# the directory to save the output of analyses
OUT_DIR=/scratch/eld72413/SAM_seq/dSNP_results

############################################
########      Dependencies GACRC     #######
############################################

# module load VEP/101.0-foss-2019b-Perl-5.30.0

module load R/4.0.0-foss-2019b

# module load BEDTools/2.30.0-GCC-8.3.0

# module load Python/3.8.6-GCCcore-10.2.0

module load BCFtools/1.13-GCC-8.3.0

module load PLINK/1.9b_5-x86_64

############################################
########      Dependencies MSI     ########
############################################

# module load R_ML/3.3.3


