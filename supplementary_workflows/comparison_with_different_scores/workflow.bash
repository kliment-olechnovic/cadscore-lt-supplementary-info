#!/bin/bash

# Set this script to fail at any command failure or any unset used variable

set -euo pipefail

# Change to the script location directory

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

################################################################################

# Main working directory

export WORKROOTDIR="$(pwd)"

# Input data directory

export INPUT_DATA_FILES_DIR="${WORKROOTDIR}/input_data"

# Software versions for download

export OLDCADSCORE_VERSION="1.1662"
export VORONOTA_VERSION="1.29.4816"
export CADSCORELT_VERSION="0.9.196"

# Software directories to create and use

export OLDCADSCORE_DIR="${WORKROOTDIR}/tools/cadscore_${OLDCADSCORE_VERSION}"
export VORONOTA_DIR="${WORKROOTDIR}/tools/voronota_${VORONOTA_VERSION}"
export CADSCORELT_DIR="${WORKROOTDIR}/tools/cadscore-lt_v${CADSCORELT_VERSION}"

# OpenStructure path to use

export OPENSTRUCTURE_APP="$HOME/miniforge3/envs/openstructure/bin/ost"

# Output data directories to create and use

export OUTPUT_DATA_CADSCORE_VARIANT_SCORES_DIR="${WORKROOTDIR}/output_data/cadscore_variant_scores"
export OUTPUT_DATA_OPENSTRUCTURE_SCORES_DIR="${WORKROOTDIR}/output_data/openstructure_scores"
export OUTPUT_DATA_OPENSTRUCTURE_AUTO_SCORES_DIR="${WORKROOTDIR}/output_data/openstructure_auto_scores"

# Output summaries directory and files to create and use

export OUTPUT_SUMMARY_DIR="${WORKROOTDIR}/output_data_summary"
export OUTPUT_SUMMARY_TABLE_OF_CADSCORE_VARIANT_SCORES="${OUTPUT_SUMMARY_DIR}/summary_table_of_cadscore_variant_scores.tsv"
export OUTPUT_SUMMARY_TABLE_OF_OPENSTRUCTURE_SCORES="${OUTPUT_SUMMARY_DIR}/summary_table_of_openstructure_scores.tsv"
export OUTPUT_SUMMARY_TABLE_OF_OPENSTRUCTURE_AUTO_SCORES="${OUTPUT_SUMMARY_DIR}/summary_table_of_openstructure_auto_scores.tsv"

# Maximum number of processors to use

export MAX_PROCESSORS_TO_USE="8"

################################################################################

# Run the workflow scripts sequentially

./scripts/get_and_check_tools.bash

./scripts/get_input_data_files.bash

./scripts/calc_cadscore_variant_scores.bash

./scripts/calc_openstructure_scores.bash

./scripts/calc_openstructure_auto_scores.bash

./scripts/analyze_summary_table_of_cadscore_variant_scores.bash

./scripts/analyze_summary_tables_of_different_scores.bash


