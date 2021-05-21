#!/usr/bin/env python3
import argparse
import gzip
import numpy as np
import pandas as pd
from libraries import codontable
import os


def open_profile(folder, ensg):
    path = "{0}/{1}_NT/sitemutsel_1.run.siteprofiles".format(folder, ensg)
    if not os.path.isfile(path): path = path.replace("_null_", "__")
    return pd.read_csv(path, sep="\t", skiprows=1, header=None,
                       names="site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y".split(","))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, nargs="+", dest="vcf",
                        help="The relative name of the .vcf file")
    parser.add_argument('--folder', required=True, type=str, dest="folder", help="The experiment folder path")
    parser.add_argument('--output', required=True, type=str, dest="output", help="The output file")
    args = parser.parse_args()

    header, profiles_dict = {}, {}
    annot_file = gzip.open(args.output, 'wt')
    duplicate = 0
    total = 0
    for vcf_path in args.vcf:
        print("Loading file " + vcf_path)
        vcf_file = gzip.open(vcf_path, 'rt')
        last_pos = -1
        last_line = ""
        for vcf_line in vcf_file:
            if vcf_line[0] == '#':
                if vcf_line[1] != '#':
                    line_list = vcf_line.strip().split("\t")
                    header = {k: i for i, k in enumerate(line_list)}
                    vcf_line = "\t".join(line_list[:header["CHR"]]) + "\n"
                annot_file.write(vcf_line)
                continue

            line_list = vcf_line.strip().split("\t")
            ensg = line_list[header["ENSG"]]
            if ensg == "None": continue

            total += 1
            if last_pos == line_list[header["POS"]]:
                print(last_line)
                print(vcf_line)
                duplicate += 1
                continue

            nuc_pos = int(line_list[header["ENSG_POS"]])
            codon_pos = int(nuc_pos / 3)

            if ensg not in profiles_dict:
                profiles_dict[ensg] = open_profile(args.folder, ensg)

            aa_ref = codontable[line_list[header["CODON_REF"]]]
            aa_alt = codontable[line_list[header["CODON_ALT"]]]
            if aa_alt == "X" or aa_ref == "X": continue
            infos = ["{0}={1}".format(k, line_list[header[k]]) for k, i in header.items() if i > header["CHR"]]
            line_list[header["INFO"]] += ";" + ";".join(infos)
            selcoeff = np.log(profiles_dict[ensg][aa_alt][codon_pos] / profiles_dict[ensg][aa_ref][codon_pos])
            line_list[header["INFO"]] += ";AA_REF={1};AA_ALT={2};SEL_COEFF={0:.3f}".format(selcoeff, aa_ref, aa_alt)
            last_line = "\t".join(line_list[:header["CHR"]]) + "\n"
            last_pos = line_list[header["POS"]]
            annot_file.write(last_line)
        vcf_file.close()

    print("{0} duplicated SNPs ({1:.2f}%)".format(duplicate, 100 * duplicate / total))
    annot_file.close()
