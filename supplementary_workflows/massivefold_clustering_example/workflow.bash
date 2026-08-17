#!/bin/bash

# Set this script to fail at any command failure or any unset used variable

set -euo pipefail

# Change to the script location directory

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

################################################################################

# Main working directory

export WORKROOTDIR="$(pwd)"

# Input data directory that already contains all the model files to cluster

export INPUT_MODEL_FILES_DIR="${WORKROOTDIR}/input_model_files/H1202"

if [ -z "$INPUT_MODEL_FILES_DIR" ] || [ ! -d "$INPUT_MODEL_FILES_DIR" ]
then
	echo >&2 "Error: the INPUT_MODEL_FILES_DIR variable in workflow.bash is not set properly"
	exit 1
fi

# Input info variables

export INPUT_INFO_DIR="${WORKROOTDIR}/input_general_info"

export INPUT_SEQUENCES_FILE="${INPUT_INFO_DIR}/sequences.fasta"
export INPUT_REFERENSE_STRUCTURE_FILE="${INPUT_INFO_DIR}/reference/8BWL-assembly1.cif"

export STOICHIOMETRY="2 2"

# Software versions for download

export VORONOTA_VERSION="1.29.4816"
export CADSCORELT_VERSION="0.9.196"

# Software directories to create and use

export CADSCORELT_DIR="${WORKROOTDIR}/tools/cadscore-lt_v${CADSCORELT_VERSION}"

# Output data directories and files to create and use

export OUTPUT_DATA_DIR="${WORKROOTDIR}/output_data"
export INPUT_MODEL_FILES_PATHS_FILE="${OUTPUT_DATA_DIR}/list_of_input_model_paths"
export OUTPUT_SUMMARY_DIR="${WORKROOTDIR}/output_data_summary"

# Maximum number of processors to use

export MAX_PROCESSORS_TO_USE="8"

################################################################################

# Run the workflow scripts sequentially

./scripts/prepare_list_of_input_model_files.bash

./scripts/get_and_check_tools.bash

./scripts/calc_all_vs_all_extremely_compact_raw.bash

./scripts/calc_all_vs_all_extremely_compact_renumbered.bash

./scripts/calc_all_vs_all_extremely_compact_renumbered_remapped_clustered.bash

./scripts/calc_reference_based_renumbered.bash

./scripts/calc_reference_based_renumbered_remapped.bash

./scripts/analyze.bash

