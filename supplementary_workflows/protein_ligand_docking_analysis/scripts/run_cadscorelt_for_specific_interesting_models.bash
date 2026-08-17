#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

OUTDIR="./output/interesting_cases"

"$CADSCORELTAPP" \
  --t ./output/prepared_input_files/3tsk/target/target_3tsk_ligand.cif \
  --m ./output/prepared_input_files/3tsk/models/model_3tsk_358.cif ./output/prepared_input_files/3tsk/models/model_3tsk_713.cif \
  --include-heteroatoms \
  --consider-residue-names \
  --processors 8 \
  --scoring-types contacts sites \
  --scoring-levels atom \
  --subselect-contacts '[-inter-chain -a1 [-chain x]]' \
  --subselect-atoms '[-protein]' \
  --output-with-identities \
  --output-with-f1 \
  --output-with-areas \
  --conflate-atom-types \
  --conflation-config-file ./output/prepared_input_files/3tsk/equivalent_atoms.txt \
  --output-dir "${OUTDIR}/results" \
  --local-output-formats table pdb contactmap graphics-pymol \
  --local-output-levels atom
