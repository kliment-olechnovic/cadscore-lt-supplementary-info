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

mapfile -t CHAINS_IDENTITY_MAP < <(awk ' /^(ATOM  |HETATM)/ {chain=substr($0, 22, 1); if (chain != " "){seen[chain] = 1}} END {for(chain in seen){print chain ":" chain}}' "$TFILE" | sort)

"$OPENSTRUCTURE_APP" compare-structures \
  -v 1 \
  -r "$TFILE" \
  -m "$MFILE" \
  -o "$OUTFILE_FULL_JSON" \
  -mf pdb \
  -rf pdb \
  -rna \
  -c "${CHAINS_IDENTITY_MAP[@]}" \
  --dockq \
  --ics \
  --ips \
  --qs-score \
  --ilddt

mkdir -p "$(dirname ${OUTFILE_SUMMARY_LOG})"

if [ ! -s "$OUTFILE_FULL_JSON" ]
then
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Failed to run openstructure for: '${TFILE}' '${MFILE}'"
	exit 1
fi

{
	echo "input ost_ilddt ost_qs_global ost_qs_best ost_ics ost_ips ost_dockq_ave" | tr ' ' '\t'
	jq -r '[.ilddt, .qs_global, .qs_best, .ics, .ips, .dockq_ave] | @tsv' "$OUTFILE_FULL_JSON" | sed "s|^|${OUTFILE_BASE}\t|"
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

