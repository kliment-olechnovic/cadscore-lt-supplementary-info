#!/bin/bash

set -u

################################################################################

SCRIPTSDIR="${WORKROOTDIR}/scripts/cadscore_variants"

OUTDIR="$1"
TFILE="$2"
MFILE="$3"

if [ -z "$OUTDIR" ] || [ ! -s "$TFILE" ] || [ ! -s "$MFILE" ]
then
	echo "Invalid input args: '${OUTDIR}' '${TFILE}' '${MFILE}'."
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
	${SCRIPTSDIR}/run_cadscore_variant_5_for_pair_of_files.bash none "${OUTDIR}/all_contacts/cadscore5" "$TFILE" "$MFILE"
	${SCRIPTSDIR}/run_cadscore_variant_5_for_pair_of_files.bash iface "${OUTDIR}/interchain_contacts/cadscore5" "$TFILE" "$MFILE"
	
	mkdir -p "$(dirname ${OUTFILE_SUMMARY_LOG})"
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Skipped fully summarizing: '${TFILE}' '${MFILE}'."
	exit 1
fi

for VERSION in 1 3 5 6 8 9
do
	${SCRIPTSDIR}/run_cadscore_variant_${VERSION}_for_pair_of_files.bash none "${OUTDIR}/all_contacts/cadscore${VERSION}" "$TFILE" "$MFILE"
done

for VERSION in 1 4 5 7 8 9
do
	${SCRIPTSDIR}/run_cadscore_variant_${VERSION}_for_pair_of_files.bash iface "${OUTDIR}/interchain_contacts/cadscore${VERSION}" "$TFILE" "$MFILE"
done

CV1_TARGET_RESIDUES="$(cat ${OUTDIR}/all_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $3}')"
CV1_TARGET_ATOMS="$(cat ${OUTDIR}/all_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV1_MODEL_RESIDUES="$(cat ${OUTDIR}/all_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $7}')"
CV1_MODEL_ATOMS="$(cat ${OUTDIR}/all_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $9}')"

CV3_TARGET_ATOMS="$(cat ${OUTDIR}/all_contacts/cadscore3/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"
CV3_MODEL_ATOMS="$(cat ${OUTDIR}/all_contacts/cadscore3/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $7}')"

CV4_TARGET_ATOMS="$(cat ${OUTDIR}/interchain_contacts/cadscore4/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $8}')"
CV4_MODEL_ATOMS="$(cat ${OUTDIR}/interchain_contacts/cadscore4/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $9}')"

CV1_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $11}')"
CV3_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore3/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $3}')"
CV5_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore5/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV5F1_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore5/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"
CV6_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore6/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $3}')"
CV8_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore8/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV8F1_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore8/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"
CV9_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore9/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV9F1_ALL_SCORE="$(cat ${OUTDIR}/all_contacts/cadscore9/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"

CV1_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore1/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $11}')"
CV4_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore4/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $3}')"
CV5_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore5/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV5F1_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore5/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"
CV7_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore7/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $3}')"
CV8_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore8/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV8F1_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore8/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"
CV9_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore9/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $5}')"
CV9F1_IFACE_SCORE="$(cat ${OUTDIR}/interchain_contacts/cadscore9/score_logs/${OUTFILE_BASE} | tail -1 | awk '{print $6}')"

CV1_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore1/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV3_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore3/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV5_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore5/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV6_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore6/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV8_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore8/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV9_ALL_TIME="$(cat ${OUTDIR}/all_contacts/cadscore9/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"

CV1_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore1/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV4_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore4/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV5_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore5/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV7_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore7/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV8_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore8/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"
CV9_IFACE_TIME="$(cat ${OUTDIR}/interchain_contacts/cadscore9/time_logs/${OUTFILE_BASE} | egrep '^\s+User time ' | awk '{print $NF}')"

readonly TMPLDIR=$(mktemp -d)
trap "rm -r $TMPLDIR" EXIT

{
cat << EOF
input
cv1_target_residues
cv1_model_residues
cv1_target_atoms
cv1_model_atoms
cv3_target_atoms
cv3_model_atoms
cv4_target_atoms
cv4_model_atoms
old_global_score
voronota_js_global_score
cadscore_lt_global_score
cadscore_lt_global_f1
cadscore_lt_rm_global_score
cadscore_lt_rm_global_f1
cadscore_lt_r_global_score
cadscore_lt_r_global_f1
old_iface_score
voronota_js_iface_score
cadscore_lt_iface_score
cadscore_lt_iface_f1
cadscore_lt_rm_iface_score
cadscore_lt_rm_iface_f1
cadscore_lt_r_iface_score
cadscore_lt_r_iface_f1
old_global_time
voronota_js_global_time
cadscore_lt_global_time
cadscore_lt_rm_global_time
cadscore_lt_r_global_time
old_iface_time
voronota_js_iface_time
cadscore_lt_iface_time
cadscore_lt_rm_iface_time
cadscore_lt_r_iface_time
voronota_js_lt_global_score
voronota_js_lt_global_time
voronota_js_lt_iface_score
voronota_js_lt_iface_time
EOF
} > "${TMPLDIR}/header"

{
cat << EOF
$OUTFILE_BASE
$CV1_TARGET_RESIDUES
$CV1_MODEL_RESIDUES
$CV1_TARGET_ATOMS
$CV1_MODEL_ATOMS
$CV3_TARGET_ATOMS
$CV3_MODEL_ATOMS
$CV4_TARGET_ATOMS
$CV4_MODEL_ATOMS
$CV1_ALL_SCORE
$CV3_ALL_SCORE
$CV5_ALL_SCORE
$CV5F1_ALL_SCORE
$CV8_ALL_SCORE
$CV8F1_ALL_SCORE
$CV9_ALL_SCORE
$CV9F1_ALL_SCORE
$CV1_IFACE_SCORE
$CV4_IFACE_SCORE
$CV5_IFACE_SCORE
$CV5F1_IFACE_SCORE
$CV8_IFACE_SCORE
$CV8F1_IFACE_SCORE
$CV9_IFACE_SCORE
$CV9F1_IFACE_SCORE
$CV1_ALL_TIME
$CV3_ALL_TIME
$CV5_ALL_TIME
$CV8_ALL_TIME
$CV9_ALL_TIME
$CV1_IFACE_TIME
$CV4_IFACE_TIME
$CV5_IFACE_TIME
$CV8_IFACE_TIME
$CV9_IFACE_TIME
$CV6_ALL_SCORE
$CV6_ALL_TIME
$CV7_IFACE_SCORE
$CV7_IFACE_TIME
EOF
} > "${TMPLDIR}/values"

mkdir -p "$(dirname ${OUTFILE_SUMMARY_LOG})"

if [ "$(cat ${TMPLDIR}/header | egrep '.' | wc -l)" != "$(cat ${TMPLDIR}/values | egrep '.' | wc -l)" ]
then
	touch "${OUTFILE_SUMMARY_LOG}"
	echo "Failed to fully summarize: '${TFILE}' '${MFILE}'."
	exit 1
fi

{
cat "${TMPLDIR}/header" | tr '\n' '\t' | sed 's/\s*$//'  | sed 's/$/\n/'
cat "${TMPLDIR}/values" | tr '\n' '\t' | sed 's/\s*$//'  | sed 's/$/\n/'
} > "${OUTFILE_SUMMARY_LOG}"

################################################################################

exit 0

