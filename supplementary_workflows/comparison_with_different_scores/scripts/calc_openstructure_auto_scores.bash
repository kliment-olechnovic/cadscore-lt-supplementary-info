#!/bin/bash

set -u

################################################################################

cat "${WORKROOTDIR}/input_data_lists/pairs_of_models_from_casp16.txt" \
| awk -v prefix="${INPUT_DATA_FILES_DIR}" '{print prefix "/" $1 " " prefix "/" $2}' \
| xargs -L 1 -P "$MAX_PROCESSORS_TO_USE" ${WORKROOTDIR}/scripts/openstructure/run_ost_auto_for_pair_of_files.bash "$OUTPUT_DATA_OPENSTRUCTURE_AUTO_SCORES_DIR"

if [ ! -s "$OUTPUT_SUMMARY_TABLE_OF_OPENSTRUCTURE_AUTO_SCORES" ]
then
	mkdir -p "$(dirname ${OUTPUT_SUMMARY_TABLE_OF_OPENSTRUCTURE_AUTO_SCORES})"
	
	find ${OUTPUT_DATA_OPENSTRUCTURE_AUTO_SCORES_DIR}/summary_logs/ -type f -not -empty \
	| sort \
	| xargs -L 1 cat \
	| awk '{if((NR==1 || $1!="input") && NF==11){print $0}}' \
	| sed 's/\s\+/\t/g' \
	> "$OUTPUT_SUMMARY_TABLE_OF_OPENSTRUCTURE_AUTO_SCORES"
fi

################################################################################

exit 0
