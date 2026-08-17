#!/bin/bash

set -u

################################################################################

MOREARGS="$1"
OUTDIR="$2"
TFILE="$3"
MFILE="$4"

if [ "$MOREARGS" == "none" ]
then
	MOREARGS=""
fi

if [ "$MOREARGS" == "iface" ]
then
	MOREARGS="--subselect-contacts [-inter-chain]"
fi

if [ -z "$OUTDIR" ] || [ ! -s "$TFILE" ] || [ ! -s "$MFILE" ]
then
	echo "Invalid input args: '${OUTDIR}' '${TFILE}' '${MFILE}' '${MOREARGS}'."
	exit 1
fi

OUTFILE_BASE="$(basename ${TFILE})__$(basename ${MFILE})"
OUTFILE_TIME_LOG=${OUTDIR}/time_logs/${OUTFILE_BASE}
OUTFILE_SCORE_LOG=${OUTDIR}/score_logs/${OUTFILE_BASE}

if [ -s "$OUTFILE_TIME_LOG" ] && [ -s "$OUTFILE_SCORE_LOG" ]
then
	exit 0
fi

readonly TMPLDIR=$(mktemp -d)
trap "rm -r $TMPLDIR" EXIT

TNAME="$(basename $(dirname ${TFILE}))"
TSEQFILE="${INPUT_DATA_FILES_DIR}/sequences/${TNAME}.fasta"

if [ -s "$TSEQFILE" ]
then
	MOREARGS="$MOREARGS --reference-sequences-file $TSEQFILE"
fi

/usr/bin/time -v \
${CADSCORELT_DIR}/cadscore-lt -t "$TFILE" -m "$MFILE" --output-with-f1 --processors 1 ${MOREARGS} \
2> "${TMPLDIR}/time_log" \
> "${TMPLDIR}/score_log"

if [ ! -s "${TMPLDIR}/time_log" ] || [ ! -s "${TMPLDIR}/score_log" ]
then
	echo "CAD-score-8 failed with: '${OUTDIR}' '${TFILE}' '${MFILE}' '${MOREARGS}'."
	exit 1
fi

mkdir -p "$(dirname ${OUTFILE_TIME_LOG})"
mv "${TMPLDIR}/time_log" "${OUTFILE_TIME_LOG}"

mkdir -p "$(dirname ${OUTFILE_SCORE_LOG})"
mv "${TMPLDIR}/score_log" "${OUTFILE_SCORE_LOG}"

################################################################################

exit 0

