#!/bin/bash

# Set this script to fail at any command failure or any unset used variable

set -euo pipefail

# Change to the script location directory

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

################################################################################

# Main working directory

export WORKROOTDIR="$(pwd)"

# Input info variables

export INPUT_INFO_DIR="${WORKROOTDIR}/input_general_info"

export INPUT_SEQUENCES_FILE="${INPUT_INFO_DIR}/sequences.fasta"

export STOICHIOMETRY="3 3 6 18 3 3"

# Input data directory

export INPUT_MODEL_FILES_DIR="${WORKROOTDIR}/input_model_files/H2424"

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

./scripts/get_and_check_tools.bash

./scripts/get_input_data_files.bash

./scripts/prepare_list_of_input_model_files.bash

./scripts/calc_all_vs_all_extremely_compact_raw.bash

./scripts/calc_all_vs_all_extremely_compact_renumbered.bash

./scripts/calc_all_vs_all_extremely_compact_renumbered_remapped_clustered.bash

./scripts/analyze.bash

