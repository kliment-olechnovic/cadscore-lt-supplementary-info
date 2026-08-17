#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

cat "./input/pdbids.txt" \
| while read -r PDBID
do
	./scripts/for_single_target/prepare_structures_for_single_target.bash "$PDBID"
	./scripts/for_single_target/run_cadscorelt_for_single_target.bash "$PDBID"
done

###########################################################################

exit 0

