#!/usr/bin/env bash
VM=tlatrill@pbil-gates.univ-lyon1.fr
DIR=/beegfs/data/tlatrill
rsync -P -v -m -r --include='*.csv' --include='*polyDFE.out' --include='*.txt' --include='*/' --exclude='*' ${VM}:${DIR}/AdaptaPop/Contrasts/ ./
rsync -P -v -m -r --include='*.pdf' --include='*.png' --include='*.tsv' --include='*/' --exclude='*' ${VM}:${DIR}/AdaptaPop/Contrasts/ ./
