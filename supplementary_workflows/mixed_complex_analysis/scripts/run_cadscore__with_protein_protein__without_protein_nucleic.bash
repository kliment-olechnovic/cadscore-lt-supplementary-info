#!/bin/bash

set -u

cd "$WORKROOTDIR"

#########################################################

OUTDIR="./output/gt_rb1_pp1_pn0_rc1e"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/targets \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
  --max-exhaustive-remappings 999 \
  --local-output-formats pdb contactmap graphics-pymol \
  --local-output-levels residue chain \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/best_models_rb1_pp1_pn0_rc1e"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/aligned/targets \
  --models ./input/aligned/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --scoring-types contacts sites \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
  --max-exhaustive-remappings 999 \
  --local-output-formats table pdb contactmap graphics-pymol \
  --local-output-levels residue chain \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/best_models_rb1_pp1_pn0_rc1"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/aligned/targets \
  --models ./input/aligned/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --scoring-types contacts sites \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
  --local-output-formats table pdb contactmap graphics-pymol \
  --local-output-levels residue chain \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_rb0_pp1_pn0_rc0"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains false \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_rb1_pp1_pn0_rc0"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains false \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_rb1_pp1_pn0_rc1"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_rb1_pp1_pn0_rc1e"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
  --max-exhaustive-remappings 999 \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_sites_rb1_pp1_pn0_rc1"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --subselect-atoms '[-protein]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-types sites \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_sites_rb1_pp1_pn0_rc1e"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --targets ./input/targets \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --subselect-atoms '[-protein]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-types sites \
  --scoring-levels residue chain \
  --output-with-f1 \
  --output-with-identities \
  --output-with-areas \
  --remap-chains true \
  --max-exhaustive-remappings 999 \
&> "${OUTDIR}/log.txt" 

#########################################################

OUTDIR="./output/run_clustering_rb1_pp1_pn0_rc1"

mkdir -p "${OUTDIR}"

[ -f "${OUTDIR}/log.txt" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --models ./input/models \
  --recursive-directory-search \
  --processors 8 \
  --reference-sequences-file ./input/sequences.fasta \
  --reference-stoichiometry 4 2 2 2 2 \
  --subselect-contacts '[-inter-chain -a1 [-protein] -a2 [-protein]]' \
  --output-dir "${OUTDIR}/results" \
  --output-global-scores _none \
  --scoring-levels residue \
  --extremely-compact-output \
  --clustering-thresholds 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 \
  --remap-chains true \
&> "${OUTDIR}/log.txt" 

#########################################################

exit 0


