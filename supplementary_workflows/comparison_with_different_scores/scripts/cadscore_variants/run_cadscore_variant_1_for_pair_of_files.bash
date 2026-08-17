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
	MOREARGS="-c"
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

/usr/bin/time -v \
${OLDCADSCORE_DIR}/bin/CADscore_calc.bash -u -b -t "$TFILE" -m "$MFILE" -D "${TMPLDIR}/db" ${MOREARGS} \
2> "${TMPLDIR}/time_log"

${OLDCADSCORE_DIR}/bin/CADscore_read_global_scores.bash -D "${TMPLDIR}/db" \
> "${TMPLDIR}/score_log"

if [ ! -s "${TMPLDIR}/time_log" ] || [ ! -s "${TMPLDIR}/score_log" ]
then
	echo "CAD-score-1 failed with: '${OUTDIR}' '${TFILE}' '${MFILE}' '${MOREARGS}'."
	exit 1
fi

mkdir -p "$(dirname ${OUTFILE_TIME_LOG})"
mv "${TMPLDIR}/time_log" "${OUTFILE_TIME_LOG}"

mkdir -p "$(dirname ${OUTFILE_SCORE_LOG})"
mv "${TMPLDIR}/score_log" "${OUTFILE_SCORE_LOG}"

################################################################################

exit 0

