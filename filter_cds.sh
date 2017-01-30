#!/usr/bin/env bash
# cd /panhome/tlatrill/AdaptaPop
cd ~/AdaptaPop
python3 ./create_interval_cds.py
bedtools intersect -a ./data/Homo_sapiens_79_polymorphism.vcf -b ./data/79_interval_cds.bed -wb > ./data/Homo_sapiens_79_polymorphism_in_cds.vcf
python3 ./mk_test_cds.py
