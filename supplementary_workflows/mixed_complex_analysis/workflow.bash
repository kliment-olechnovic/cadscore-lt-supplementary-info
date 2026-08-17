#!/bin/bash

# Set this script to fail at any command failure or any unset used variable

set -euo pipefail

# Change to the script location directory

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

#########################################################

# Main working directory

export WORKROOTDIR="$(pwd)"

# Software versions for download

export VORONOTA_VERSION="1.29.4816"
export CADSCORELT_VERSION="0.9.196"

# Software directories to create and use

export CADSCORELT_DIR="${WORKROOTDIR}/tools/cadscore-lt_v${CADSCORELT_VERSION}"

export CADSCORELTAPP="./tools/cadscore-lt_v${CADSCORELT_VERSION}/cadscore-lt"

#########################################################

# Run the workflow scripts sequentially

./scripts/get_input_data_files.bash

./scripts/get_and_check_tools.bash

./scripts/run_cadscore__without_protein_protein__with_protein_nucleic.bash

./scripts/run_cadscore__with_protein_protein__without_protein_nucleic.bash

./scripts/run_cadscore__with_protein_protein__with_protein_nucleic.bash

./scripts/analyze.bash

