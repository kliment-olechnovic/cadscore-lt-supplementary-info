#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

if [ ! -d "./workvenv" ]
then
	python3 -m venv workvenv

	source ./workvenv/bin/activate

	pip install rdkit

	pip install gemmi
fi

