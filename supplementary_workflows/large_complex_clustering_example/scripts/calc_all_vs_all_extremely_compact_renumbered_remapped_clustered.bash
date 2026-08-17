#!/bin/bash

set -u

################################################################################

OUTDIR="${OUTPUT_DATA_DIR}/extremely_compact/renumbered_remapped_clustered"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"${CADSCORELT_DIR}/cadscore-lt" \
  --processors "$MAX_PROCESSORS_TO_USE" \
  --reference-stoichiometry $STOICHIOMETRY \
  --reference-sequences-file "$INPUT_SEQUENCES_FILE" \
  --subselect-contacts '[-inter-chain]' \
  --consider-residue-names \
  --remap-chains \
  --extremely-compact-output \
  --clustering-thresholds 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 \
  --output-dir "${OUTDIR}/results" \
< "$INPUT_MODEL_FILES_PATHS_FILE" \
&> "${OUTDIR}/log.txt"

################################################################################

exit 0
