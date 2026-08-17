#!/bin/bash

set -u

################################################################################

OUTDIR="$1"
TFILE="$2"
MFILE="$3"

if [ -z "$OUTDIR" ] || [ ! -s "$TFILE" ] || [ ! -s "$MFILE" ]
then
	echo "Invalid input args: '${OUTDIR}' '${TFILE}' '${MFILE}'"
	exit 1
fi

OUTFILE_BASE="$(basename ${TFILE})__$(basename ${MFILE})"
OUTFILE_SUMMARY_LOG=${OUTDIR}/summary_logs/${OUTFILE_BASE}

if [ -f "$OUTFILE_SUMMARY_LOG" ]
then
	exit 0
fi

TOTAL_ATOMS="$(cat ${TFILE} ${MFILE} | egrep '^ATOM' | wc -l)"

if [ "$TOTAL_ATOMS" -gt "100000" ]
then
	mkdir -p "$(dirname ${OUTFILE_SUMMARY_LOG})"
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Skipped fully summarizing: '${TFILE}' '${MFILE}'"
	exit 1
fi

OUTFILE_FULL_JSON=${OUTDIR}/full_results/${OUTFILE_BASE}.json
mkdir -p "$(dirname ${OUTFILE_FULL_JSON})"

"$OPENSTRUCTURE_APP" compare-structures \
  -v 1 \
  -r "$TFILE" \
  -m "$MFILE" \
  -o "$OUTFILE_FULL_JSON" \
  -mf pdb \
  -rf pdb \
  -rna \
  --dockq \
  --ics \
  --ips \
  --qs-score \
  --ilddt \
  --tm-score

mkdir -p "$(dirname ${OUTFILE_SUMMARY_LOG})"

if [ ! -s "$OUTFILE_FULL_JSON" ]
then
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Failed to run openstructure for: '${TFILE}' '${MFILE}'"
	exit 1
fi

{
	echo "input ost_ilddt ost_qs_global ost_qs_best ost_ics ost_ips ost_dockq_ave ost_dockq_wave ost_dockq_ave_full ost_dockq_wave_full ost_tm_score" | tr ' ' '\t'
	jq -r '[.ilddt, .qs_global, .qs_best, .ics, .ips, .dockq_ave, .dockq_wave, .dockq_ave_full, .dockq_wave_full, .tm_score] | @tsv' "$OUTFILE_FULL_JSON" | sed "s|^|${OUTFILE_BASE}\t|"
} \
> "$OUTFILE_SUMMARY_LOG"

if [ "$(cat ${OUTFILE_SUMMARY_LOG} | awk '{print $2}' | grep -E '\S' | wc -l)" -ne "2" ]
then
	rm -f "${OUTFILE_SUMMARY_LOG}"
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Failed to collect openstructure scores for: '${TFILE}' '${MFILE}'"
	exit 1
fi

################################################################################

exit 0

