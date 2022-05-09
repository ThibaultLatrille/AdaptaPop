#!/usr/bin/env python3
import argparse
import gzip
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple
from scipy.optimize import brentq
from libraries import codontable, snp_data_frame, build_divergence_dico
import os


def s_to_omega(s):
    if s < -10:
        return -s * np.exp(s)
    else:
        return s / (1 - np.exp(-s))


def omega_to_s(omega):
    if omega == 0.0:
        return -np.inf
    a, b = -100, 1.0
    if omega > 1:
        b = omega + 1.0
    if s_to_omega(a) > omega:
        print("a", a, s_to_omega(a), omega)
    elif s_to_omega(b) < omega:
        print("b", b, s_to_omega(b), omega)
    s = brentq(lambda x: s_to_omega(x) - omega, a=a, b=b)
    return s


def open_profile(folder, ensg):
    path = "{0}/{1}_NT/sitemutsel_1.run.siteprofiles".format(folder, ensg)
    if not os.path.isfile(path):
        path = path.replace("_null_", "__")
    return pd.read_csv(path, sep="\t", skiprows=1, header=None,
                       names="site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y".split(","))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="The relative name of the .vcf file")
    parser.add_argument('--tsv', required=True, type=str, dest="tsv", help="The relative name of the .tsv file")
    parser.add_argument('--folder', required=True, type=str, dest="folder", help="The experiment folder path")
    parser.add_argument('--output', required=True, type=str, dest="output", help="The output file")
    args = parser.parse_args()
    SNP_row = namedtuple('SNP_row', ['p', 'anc', 'der', 'polarized', 'count', 'count_polarized', 'sample_size'])

    dico_snps = defaultdict(dict)
    df_snp, _, _ = snp_data_frame(args.tsv, polarize_snps=True, remove_fixed=False)
    for ensg, ddf in df_snp:
        for _, row in ddf.iterrows():
            dico_snps[ensg][row["NUC_POS"]] = SNP_row(row["ANC_PROBA"], row["ANC"], row["DER"], row["POLARIZED"],
                                                      row["COUNT"], row["COUNT_POLARIZED"], row["SAMPLE_SIZE"])

    header, profiles_dict = {}, {}
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, list(df_snp.groups), gene_level=False)
    annot_file = gzip.open(args.output, 'wt')
    duplicate = 0
    total = 0
    print("Loading file " + args.vcf)
    vcf_file = gzip.open(args.vcf, 'rt')
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

        if "CODON_REF" not in header:
            break

        line_list = vcf_line.strip().split("\t")
        ensg = line_list[header["ENSG"]]
        if ensg == "None":
            continue

        total += 1
        if last_pos == line_list[header["POS"]]:
            duplicate += 1
            continue

        nuc_pos = int(line_list[header["ENSG_POS"]])
        codon_pos = int(nuc_pos / 3)

        if nuc_pos not in dico_snps[ensg]:
            continue

        if ensg not in profiles_dict:
            profiles_dict[ensg] = open_profile(args.folder, ensg)

        codon_ref, codon_alt = line_list[header["CODON_REF"]], line_list[header["CODON_ALT"]]
        aa_ref, aa_alt = codontable[codon_ref], codontable[codon_alt]
        if aa_alt == "X" or aa_ref == "X":
            continue

        if dico_snps[ensg][nuc_pos].polarized:
            aa_anc, aa_der = aa_ref, aa_alt
            codon_anc, codon_der = codon_ref, codon_alt
            assert dico_snps[ensg][nuc_pos].anc == line_list[header["ENSG_REF"]]
        else:
            aa_anc, aa_der = aa_alt, aa_ref
            codon_anc, codon_der = codon_alt, codon_ref

        selcoeff = np.log(profiles_dict[ensg][aa_der][codon_pos] / profiles_dict[ensg][aa_anc][codon_pos])
        infos = ["{0}={1}".format(k, line_list[header[k]]) for k, i in header.items() if i > header["CHR"]]
        last_pos = line_list[header["POS"]]

        line_list[header["INFO"]] += f";ANC_PROBA={dico_snps[ensg][nuc_pos].p}"
        line_list[header["INFO"]] += f";SAMPLE_SIZE={dico_snps[ensg][nuc_pos].sample_size}"
        line_list[header["INFO"]] += f";NUC_ANC={dico_snps[ensg][nuc_pos].anc}"
        line_list[header["INFO"]] += f";NUC_DER={dico_snps[ensg][nuc_pos].der}"
        line_list[header["INFO"]] += f";CODON_ANC={codon_anc};CODON_DER={codon_der}"
        line_list[header["INFO"]] += f";AA_REF={aa_ref};AA_ALT={aa_alt}"
        line_list[header["INFO"]] += f";AA_ANC={aa_anc};AA_DER={aa_der};SEL_COEFF={selcoeff:.3f}"
        line_list[header["INFO"]] += f";POLARIZED={dico_snps[ensg][nuc_pos].polarized}"
        line_list[header["INFO"]] += f";COUNT={dico_snps[ensg][nuc_pos].count}"
        line_list[header["INFO"]] += f";COUNT_POLARIZED={dico_snps[ensg][nuc_pos].count_polarized}"
        line_list[header["INFO"]] += f";SAMPLE_SIZE={dico_snps[ensg][nuc_pos].sample_size}"
        line_list[header["INFO"]] += ";" + ";".join(infos)

        site_omega = dico_omega[ensg][codon_pos][1]
        line_list[header["INFO"]] += ";SITE_OMEGA=" + str(site_omega)
        line_list[header["INFO"]] += ";SITE_OMEGA_SEL_COEFF=" + str(omega_to_s(site_omega))

        site_omega_0 = dico_omega_0[ensg][codon_pos][1]
        line_list[header["INFO"]] += ";SITE_OMEGA_0=" + str(site_omega_0)
        line_list[header["INFO"]] += ";SITE_OMEGA_0_SEL_COEFF=" + str(omega_to_s(site_omega_0))
        last_line = "\t".join(line_list[:header["CHR"]]) + "\n"
        annot_file.write(last_line)
    vcf_file.close()

    if total != 0:
        print("{0} duplicated SNPs ({1:.2f}%)".format(duplicate, 100 * duplicate / total))
    annot_file.close()
