#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=6gb
#SBATCH --tmp=2gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edittmar@umn.edu
#SBATCH -p small,ram256g,ram1t,max
#SBATCH -o bad_mut_compile_predictions.sh.%j.out
#SBATCH -e bad_mut_compile_predictions.sh.%j.err

### *** NOTE: This script was written by Chaochih Liu for Barley project ****

set -e
set -o pipefail

# This script compiles all predictions in parallel into a single file for easy downstream processing.
#   BAD_Mutations.py will output the compiled report in the predict output subdirectories. These will then be compiled into a single report in ${OUT_DIR}.

# Usage: sbatch bad_mut_compile_predict-${PROJECT}.job

# Dependencies
module load parallel
module load python3/3.6.3_anaconda5.0.1

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User provided input arguments
# Full path to a list of predict output directories that contain *_Predictions.txt files
#   This can be generated as follows from within the predict script ${OUT_DIR}:
#   cd /path/to/predict_output_${SAMPLE_NAME}
#   Generate the list of output directories
#   find $(pwd -P) -maxdepth 1 -type d -name "*list*" > predict_out_dir_list.txt
#   Get path to list of output directories
#   find $(pwd -P) -name predict_out_dir_list.txt
PREDICT_DIR_LIST=/scratch.global/edittmar/predict_output_Sunflower_SAM/predict_out_dir_list.txt

# Full path to the list of long substitutions .txt file
#   i.e., The long_substitutions.txt file output from the script VeP_to_Subs.py
LONG_SUBS_FILE=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/SAM_SNP_BadMut_Summary_edit
# Project name will be used as prefix for final combined report
PROJECT=Sunflower_SAM

# Where do we want to store our final compiled predictions report?
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/CompiledResults

# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT=/panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations/BAD_Mutations.py

# Full path to the script bad_mut_compile_predictions.sh
COMPILE_SCRIPT=/panfs/roc/groups/9/morrellp/edittmar/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/bad_mut_compile_predictions.sh

#------------------------------
# Run compile predictions script
${COMPILE_SCRIPT} \
    ${PREDICT_DIR_LIST} \
    ${LONG_SUBS_FILE} \
    ${PROJECT} \
    ${OUT_DIR} \
    ${BAD_MUT_SCRIPT}
