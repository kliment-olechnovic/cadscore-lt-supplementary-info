#!/bin/bash

set -u

################################################################################

OUTDIR="${OUTPUT_DATA_DIR}/extremely_compact/renumbered"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"${CADSCORELT_DIR}/cadscore-lt" \
  --processors "$MAX_PROCESSORS_TO_USE" \
  --reference-stoichiometry $STOICHIOMETRY \
  --reference-sequences-file "$INPUT_SEQUENCES_FILE" \
  --subselect-contacts '[-inter-chain]' \
  --consider-residue-names \
  --extremely-compact-output \
  --output-dir "${OUTDIR}/results" \
< "$INPUT_MODEL_FILES_PATHS_FILE" \
&> "${OUTDIR}/log.txt"

################################################################################

exit 0
