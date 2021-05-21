#!/usr/bin/env bash
#### Local
VM=ubuntu@X.X.X.X
ssh -X ${VM}

#### Remote
cd ~/data/mydatalocal/
git clone https://github.com/ThibaultLatrille/AdaptaPop
cd AdaptaPop
mkdir ./OrthoMam/Datasets
chmod 755 install.sh
./install.sh

#### Local (copy dataset to remote)
cd ~/AdaptaPop
scp Polymorphism/CDS.ANNOT.tar.xz ${VM}:~/data/mydatalocal/AdaptaPop/Polymorphism
scp OrthoMam/Experiments.tar.xz ${VM}:~/data/mydatalocal/AdaptaPop/OrthoMam
scp OrthoMam/Datasets/omm.tar.xz ${VM}:~/data/mydatalocal/AdaptaPop/OrthoMam/Datasets

#### Remote (extract datasets)
cd Polymorphism
tar -xf CDS.ANNOT.tar.xz
cd ../OrthoMam
tar -xf Experiments.tar.xz
cd Datasets
tar -xf omm.tar.xz
screen -dmS Snakemake bash -c "cd ~/data/mydatalocal/AdaptaPop/Contrasts && snakemake -j 32 -k --printshellcmds"

#### Local (copy results from remote)
scp -r ${VM}:~/data/mydatalocal/AdaptaPop/Contrasts Contrasts
rsync ${VM}:~/data/mydatalocal/AdaptaPop/Contrasts/*/histogram.pdf
scp ${VM}:~/data/mydatalocal/AdaptaPop/Contrasts/*/histogram.pdf ./Contrasts/Chlorocebus_sabaeus-Barbados-gene/histogram.pdf