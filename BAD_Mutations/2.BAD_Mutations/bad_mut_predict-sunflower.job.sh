#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10gb
#SBATCH --tmp=2gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edittmar@umn.edu
#SBATCH -p small,ram256g,ram1t,max
#SBATCH -o bad_mut_predict.sh.%A_%a.out
#SBATCH -e bad_mut_predict.sh.%A_%a.err

### *** NOTE: This script was written by Chaochih Liu for Barley project ****

set -e
set -o pipefail

# This script stores all filepaths and calls on the script bad_mut_predict.sh

# Dependencies
module load parallel
module load python3/3.6.3_anaconda5.0.1

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User provided input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
# List of lists here should only include FASTA files that had an alignment
#	(i.e., .fa and .tree files were generated) in the previous align step
FASTA_LIST_OF_LISTS=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/fasta_lists/all_cds_Hannuus_list_of_lists.txt

# Full path to the config file
CONFIG_FILE=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/config.txt

# Full path to a list of MSA_Output directories that contain *.fa and *.tree files
MSA_DIR_LIST=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/MSA_output/MSA_output_dir_list.txt

# Full path to per transcript substitutions .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/BadMutationsSubs
#/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/SAM_SNP_BadMut_Summary_edit

# Full path to a list of primary transcripts, one per line
PRIMARY_TRANSCRIPTS=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Transcript_names.txt

# Sample name will be used as a prefix for outputs
SAMPLE_NAME="Sunflower_SAM"

# Full path to output directory
OUT_DIR=/scratch.global/edittmar/predict_output_${SAMPLE_NAME}
#liux1299/bad_mutations/predict_output_${SAMPLE_NAME}

# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT=/panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations/BAD_Mutations.py

# Full path to where we want to store the log files output from parallel
LOG_FILE_DIR=${OUT_DIR}/all_log_files

# Script that stores predict_sub function
PREDICT_SCRIPT=/panfs/roc/groups/9/morrellp/edittmar/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/bad_mut_predict.sh

#------------------------------
# Run predict script
${PREDICT_SCRIPT} ${FASTA_LIST_OF_LISTS} \
	${CONFIG_FILE} \
	${MSA_DIR_LIST} \
	${SUBS_DIR} \
	${PRIMARY_TRANSCRIPTS} \
	${SAMPLE_NAME} \
	${OUT_DIR} \
	${BAD_MUT_SCRIPT} \
	${LOG_FILE_DIR}
