cd /panhome/tlatrill
python ./create_interval_cds.py
bedtools intersect -a ./Homo_sapiens_79_polymorphism.vcf -b ./79_interval_cds.bed -wb > ./Homo_sapiens_79_polymorphism_in_cds.vcf
