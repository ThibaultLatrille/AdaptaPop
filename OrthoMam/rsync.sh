#!/usr/bin/env bash
VM=tlatrill@pbil-gates.univ-lyon1.fr
DIR=/beegfs/data/tlatrill
rsync -r --include='*.run.ci0.*' --include='*.run.siteprofiles'  --include='*.ali'  --include='*.rootree' --include='*/' --exclude='*' ${VM}:${DIR}/AdaptaPop/OrthoMam/Experiments/ Experiments/
