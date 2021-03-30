#!/usr/bin/env bash

VCF="AUCH.genus_snps.CHIR1_0.20140928"
VCF_TARGET="AUCH.genus_snps.ARS1"
TARGET='GCA_000317765.1_CHIR_1.0_genomic.fna'

CrossMap.py vcf flo/combined.chn.sorted ${VCF}.vcf.gz flo/"${TARGET}" "${VCF_TARGET}".vcf.gz
