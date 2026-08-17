#!/bin/bash

set -u

################################################################################

mkdir -p "$(dirname ${INPUT_MODEL_FILES_PATHS_FILE})"

if [ ! -s "$INPUT_MODEL_FILES_PATHS_FILE" ]
then
	{
		find "$INPUT_MODEL_FILES_DIR" -type f -name '*.pdb'
		find "$INPUT_MODEL_FILES_DIR" -type f -name '*.cif'
	} \
	> "$INPUT_MODEL_FILES_PATHS_FILE"
fi

if [ ! -s "$INPUT_MODEL_FILES_PATHS_FILE" ]
then
	echo >&2 "Error: no .pdb or .cif file found in the provided input directory INPUT_MODEL_FILES_DIR=${INPUT_MODEL_FILES_DIR}"
	exit 1
fi

################################################################################

exit 0
