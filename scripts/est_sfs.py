#!/usr/bin/env python3
import argparse
import gzip
import os
import pandas as pd
from libraries import complement, nucleotides

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=True, type=str, dest="vcf", help="The relative name of the .vcf file")
    parser.add_argument('-e', '--exec', required=True, type=str, dest="exec", help="The executable path")
    parser.add_argument('-c', '--config', required=True, type=str, dest="config", help="The config path")
    parser.add_argument('-s', '--seed', required=True, type=str, dest="seed", help="The seed path")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="The output file")
    args = parser.parse_args()

    print("Loading file " + args.vcf)
    fixed_poly = dict()
    snp_table, header = {"CHR": [], "REF": [], "ALT": [], "ANC": [], "COUNT": [], "SAMPLE_SIZE": [],
                         "ENSG": [], "POS": [], "TYPE": []}, {}
    vcf_file = gzip.open(args.vcf, 'rt')
    est_sfs_path = args.output.replace(".tsv.gz", ".est-sfs")
    est_sfs_file = open(est_sfs_path + ".txt", 'w')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                header = {k: i for i, k in enumerate(vcf_line.strip().split("\t"))}
            continue

        line_list = vcf_line.strip().split("\t")
        ensg = line_list[header["ENSG"]]
        if ensg == "None": continue

        if line_list[header["CHR"]] in ["X", "Y", "MT"]: continue

        genotypes = [s.split(":")[0].count("1") for s in line_list[header["FORMAT"] + 1:header["CHR"]] if
                     ("|" in s) or ("/" in s)]
        count = sum(genotypes)
        n = len(genotypes) * 2
        if count == 0:
            continue

        nuc_pos = int(line_list[header["ENSG_POS"]])
        codon_pos = int(nuc_pos / 3)
        nuc_ref = line_list[header["REF"]]
        nuc_alt = line_list[header["ALT"]]
        if line_list[header["STRAND"]] == "-":
            nuc_ref = complement[nuc_ref]
            nuc_alt = complement[nuc_alt]

        assert nuc_ref == line_list[header["ENSG_REF"]]
        sfs = [list(), list(), list(), list()]
        for nuc in nucleotides:
            if nuc == nuc_alt:
                sfs[0].append(count)
            elif nuc == nuc_ref:
                sfs[0].append(n - count)
            else:
                sfs[0].append(0)
        for gr in [1, 2, 3]:
            out_nuc = line_list[header["OUTGROUP_{0}".format(gr)]]
            for nuc in nucleotides:
                if nuc == out_nuc:
                    sfs[gr].append(1)
                else:
                    sfs[gr].append(0)

        est_sfs_file.write("\t".join([",".join(map(str, s)) for s in sfs]) + "\n")
        snp_table["CHR"].append(line_list[header["CHR"]])
        snp_table["REF"].append(nuc_ref)
        snp_table["ALT"].append(nuc_alt)
        snp_table["COUNT"].append(count)
        snp_table["SAMPLE_SIZE"].append(n)
        snp_table["ENSG"].append(ensg)
        snp_table["POS"].append(codon_pos)
        snp_table["TYPE"].append(line_list[header["SNP_TYPE"]])

    est_sfs_file.close()
    vcf_file.close()
    print(args.vcf + " loaded.")
    os.system("{0} {1} {2}.txt {3} {2}.usfs {2}.pvalues.txt".format(args.exec, args.config, est_sfs_path, args.seed))

    sfs_outfile = open(est_sfs_path + ".pvalues.txt", 'r')
    for sfs_line in sfs_outfile:
        line_split = sfs_line.strip().split(" ")
        if line_split[0] == "0": continue
        anc = ""
        max_p = 0.0
        for nuc, p_str in zip(nucleotides, line_split[3:]):
            p = float(p_str)
            if p >= max_p:
                anc = nuc
                max_p = p
        snp_table["ANC"].append(anc)
    pd.DataFrame(snp_table).to_csv(args.output, index=False, compression="gzip")
