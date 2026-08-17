#!/bin/bash

set -u

################################################################################

OUTDIR="${OUTPUT_DATA_DIR}/extremely_compact/raw"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"${CADSCORELT_DIR}/cadscore-lt" \
  --processors "$MAX_PROCESSORS_TO_USE" \
  --subselect-contacts '[-inter-chain]' \
  --extremely-compact-output \
  --output-dir "${OUTDIR}/results" \
< "$INPUT_MODEL_FILES_PATHS_FILE" \
&> "${OUTDIR}/log.txt"

################################################################################

exit 0
