#!/usr/bin/env bash
#### Local
REMOTE=ubuntu@X.X.X.X
REMOTE_FOLDER=/scratch/tlatrill/
LOCAL_FOLDER=~/Documents
ssh -X ${REMOTE}

#### Remote
cd ${REMOTE_FOLDER}
git clone https://github.com/ThibaultLatrille/AdaptaPop
cd ${REMOTE_FOLDER}/AdaptaPop
mkdir ./Contrasts/pickle
mkdir ./OrthoMam/Datasets
chmod 755 install.sh
./install.sh

#### Local (copy dataset to remote)
cd ${LOCAL_FOLDER}/AdaptaPop
scp utils.tar.xz ${REMOTE}:${REMOTE_FOLDER}/AdaptaPop
scp Polymorphism/CDS.ANNOT.tar.xz ${REMOTE}:${REMOTE_FOLDER}/AdaptaPop/Polymorphism
scp OrthoMam/Experiments.tar.xz ${REMOTE}:${REMOTE_FOLDER}/AdaptaPop/OrthoMam
scp OrthoMam/Datasets/omm.tar.xz ${REMOTE}:${REMOTE_FOLDER}/AdaptaPop/OrthoMam/Datasets

#### Remote (extract datasets)
cd ${REMOTE_FOLDER}/AdaptaPop
tar -xf utils.tar.xz
cd ${REMOTE_FOLDER}/AdaptaPop/Polymorphism
tar -xf CDS.ANNOT.tar.xz
cd ${REMOTE_FOLDER}/AdaptaPop/OrthoMam
tar -xf Experiments.tar.xz
cd ${REMOTE_FOLDER}/AdaptaPop/OrthoMam/Datasets
tar -xf omm.tar.xz
screen -dmS Snakemake bash -c "cd ${REMOTE_FOLDER}/AdaptaPop/Contrasts && snakemake -j 32 -k --printshellcmds"


#### Compress and copy experimetn
EXP="polyDFE_modelC_no_control"
#### Remote (compress results)
cd ${REMOTE_FOLDER}/AdaptaPop/Contrasts
rm -rf ./${EXP}/*/tmp/
tar -zcvf "${EXP}.tar.gz" "${EXP}"

#### Local (copy results from remote)
scp -r ${REMOTE}:${REMOTE_FOLDER}/AdaptaPop/Contrasts/${EXP} ${LOCAL_FOLDER}/AdaptaPop/Contrasts
cd ${LOCAL_FOLDER}/AdaptaPop/Contrasts
gunzip "${EXP}.tar.gz"