#!/usr/bin/env bash
git clone https://github.com/ThibaultLatrille/AdaptaPop
cd AdaptaPop
sudo apt install -qq -y python3-dev python3-pip
pip3 install snakemake scipy numpy matplotlib pandas ete3 bio statsmodels --user
sudo apt install -qq -y bedtools paml
sudo apt install -qq -y make cmake clang
cd OrthoMam
git clone https://github.com/ThibaultLatrille/bayescode && cd bayescode && git checkout dev && make release && cd ../..
