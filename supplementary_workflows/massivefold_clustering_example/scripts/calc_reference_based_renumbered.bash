#!/bin/bash

set -u

################################################################################

OUTDIR="${OUTPUT_DATA_DIR}/normal/reference_based_renumbered"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"${CADSCORELT_DIR}/cadscore-lt" \
  -t "$INPUT_REFERENSE_STRUCTURE_FILE" \
  --processors "$MAX_PROCESSORS_TO_USE" \
  --reference-stoichiometry $STOICHIOMETRY \
  --reference-sequences-file "$INPUT_SEQUENCES_FILE" \
  --subselect-contacts '[-inter-chain]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
< "$INPUT_MODEL_FILES_PATHS_FILE" \
&> "${OUTDIR}/log.txt"

################################################################################

exit 0
