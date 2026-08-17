#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

PDBID="$1"

ROOT_INDIR="./output/prepared_input_files/${PDBID}"
TARGET_INDIR="${ROOT_INDIR}/target"
MODELS_INDIR="${ROOT_INDIR}/models"

TARGET_INFILE="$(find ${TARGET_INDIR} -type f)"

if [ ! -s "$TARGET_INFILE" ]
then
	echo "No target file ${TARGET_INFILE}."
	exit 1
fi

if [ ! -d "$MODELS_INDIR" ]
then
	echo "No models directory ${MODELS_INDIR}."
	exit 1
fi

OUTDIR="./output/cadscore/${PDBID}"

if [ -d "$OUTDIR" ]
then
	echo "Skipping already assessed ${PDBID}."
	exit 0
fi

mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/global_scores_raw.tsv"
LOGFILE="${OUTDIR}/log_raw.txt"

[ -f "$OUTFILE" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --t "$TARGET_INFILE" \
  --m "$MODELS_INDIR" \
  --include-heteroatoms \
  --consider-residue-names \
  --processors 8 \
  --scoring-types contacts sites \
  --scoring-levels atom residue \
  --subselect-contacts '[-inter-chain -a1 [-chain x]]' \
  --subselect-atoms '[-protein]' \
  --output-with-identities \
  --output-with-f1 \
  --output-with-areas \
  --output-global-scores "$OUTFILE" \
&> "$LOGFILE"

LEQFILE="${ROOT_INDIR}/equivalent_atoms.txt"

CONFLATION_OPTIONS=" --conflate-atom-types --conflation-config-file $LEQFILE "

if [ ! -s "$LEQFILE" ]
then
	CONFLATION_OPTIONS=" --conflate-atom-types "
fi

OUTFILE="${OUTDIR}/global_scores_with_atom_types_info.tsv"
LOGFILE="${OUTDIR}/log_with_atom_types_info.txt"

[ -f "$OUTFILE" ] || \
/usr/bin/time -v \
"$CADSCORELTAPP" \
  --t "$TARGET_INFILE" \
  --m "$MODELS_INDIR" \
  --include-heteroatoms \
  --consider-residue-names \
  --processors 8 \
  --scoring-types contacts sites \
  --scoring-levels atom residue \
  --subselect-contacts '[-inter-chain -a1 [-chain x]]' \
  --subselect-atoms '[-protein]' \
  --output-with-identities \
  --output-with-f1 \
  --output-with-areas \
  $CONFLATION_OPTIONS \
  --output-global-scores "$OUTFILE" \
&> "$LOGFILE"

