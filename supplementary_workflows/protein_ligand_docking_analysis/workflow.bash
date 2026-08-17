#!/bin/bash

# Set this script to fail at any command failure or any unset used variable

set -euo pipefail

# Change to the script location directory

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

#########################################################

# Main working directory

export WORKROOTDIR="$(pwd)"

# CASF-2016 input directory

export CASF2016_INPUT_DIR="/data/datasets/CASF16/CASF-2016"

if [ -z "$CASF2016_INPUT_DIR" ] || [ ! -d "$CASF2016_INPUT_DIR" ]
then
	echo >&2 "Error: the CASF2016_INPUT_DIR variable in workflow.bash is not set properly"
	exit 1
fi

# Software versions for download

export VORONOTA_VERSION="1.29.4816"
export CADSCORELT_VERSION="0.9.196"
export CADSCORELT_EXTRAUTILS_VERSION="0.9.1"

# Software directories to create and use

export CADSCORELT_DIR="${WORKROOTDIR}/tools/cadscore-lt_v${CADSCORELT_VERSION}"
export CADSCORELTAPP="${WORKROOTDIR}/tools/cadscore-lt_v${CADSCORELT_VERSION}/cadscore-lt"

export CADSCORELT_EXTRAUTILS_DIR="${WORKROOTDIR}/tools/extras/cadscorelt-extra-data-utilities-${CADSCORELT_EXTRAUTILS_VERSION}"

#########################################################

# Run the workflow scripts sequentially

./scripts/get_and_check_tools.bash

./scripts/ensure_venv.bash

./scripts/process_all_input_target_pdb_ids.bash

./scripts/run_cadscorelt_for_specific_interesting_models.bash

./scripts/analyze.bash

