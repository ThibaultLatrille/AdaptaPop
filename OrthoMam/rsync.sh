#!/usr/bin/env bash
rsync -r --include='*.run.ci0.*' --exclude="*.chain.gz" tlatrill@pbil-gates.univ-lyon1.fr:/beegfs/data/tlatrill/AdaptaPop/OrthoMam/Experiments/ Experiments/
