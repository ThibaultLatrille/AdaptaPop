#!/usr/bin/env bash
# cd /panhome/tlatrill/AdaptaPop
cd ~/AdaptaPop
python3 ./write_bedfile.py
bedtools intersect -a ./data/Homo_sapiens_88_GRCh38_polymorphism.vcf -b ./data/88_GRCh38_interval_cds.bed -wb > ./data/Homo_sapiens_88_GRCh38_polymorphism_in_cds.vcf

python3 ./write_estimates_snp.py
