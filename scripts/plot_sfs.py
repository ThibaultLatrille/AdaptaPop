import argparse
import gzip
import numpy as np
from libraries import RED, GREEN, BLUE, LIGHTGREEN, YELLOW
from collections import defaultdict
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 256
weak = ["A", "T"]
strong = ["G", "C"]
nbr_replicates = 5
subsample = 16


def daf_to_count(daf_list, min_n):
    array_daf = np.array(daf_list, dtype=np.int64).T
    return np.array([np.array(np.bincount(d, minlength=min_n)[1:]) / sum([1 for i in d if i != 0]) for d in array_daf])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=False, type=str, dest="vcf", help="Input vcf file")
    parser.add_argument('-s', '--method', required=False, type=str, dest="method", help="Sel coeff parameter")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output pdf file")
    args = parser.parse_args()

    assert args.method in ["MutSel", "Omega", "Omega_0", "WS", "SW", "SIFT"]

    vcf_file = gzip.open(args.vcf, 'rt')
    list_non_syn, list_syn_daf = [], []
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                line_strip = vcf_line.strip()
                if args.method == "SIFT":
                    line_strip = line_strip.replace("\t", "\tSIFT_")
                    line_strip += "\t#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                header = {k: i for i, k in enumerate(line_strip.split("\t"))}
            continue

        split_line = vcf_line.strip().split("\t")
        info = split_line[header["INFO"]]
        dico_info = {k: v for k, v in [s.split("=") for s in info.split(";") if "=" in s]}

        if args.method == "WS" and not ((dico_info["ANC"] in weak) and (dico_info["DER"] in strong)):
            continue

        if args.method == "SW" and not ((dico_info["ANC"] in strong) and (dico_info["DER"] in weak)):
            continue

        polarized = dico_info["POLARIZED"] == "True"
        if not polarized and args.method == "SIFT":
            continue

        n = int(dico_info["SAMPLE_SIZE"])
        k = int(dico_info["COUNT_POLARIZED"])
        max_daf = min(n, subsample)
        if n > subsample:
            count = [np.random.hypergeometric(k, n - k, subsample) for i in range(nbr_replicates)]
            count = [i if i != 0 and i != max_daf else 0 for i in count]
        else:
            count = [k]

        if np.all(np.isnan(count)):
            continue

        if dico_info["SNP_TYPE"] == "Syn":
            list_syn_daf.append(count)
        elif dico_info["SNP_TYPE"] == "NonSyn":
            if args.method in ["MutSel", "WS", "SW", "SIFT"]:
                S = float(dico_info["SEL_COEFF"])
            elif args.method == "Omega":
                S = float(dico_info["SITE_OMEGA_SEL_COEFF"])
            else:
                assert args.method == "Omega_0"
                S = float(dico_info["SITE_OMEGA_0_SEL_COEFF"])

            if args.method == "SIFT":
                sif_info = split_line[header["SIFT_INFO"]]
                sift_find = sif_info.find("Sift=0")
                if sift_find == -1:
                    continue
                sift_score = float(sif_info[sift_find:].split("|")[2])
                list_non_syn.append((sift_score, S, count))
            elif np.isfinite(S):
                list_non_syn.append((S, count))

    vcf_file.close()

    syn_count = daf_to_count(list_syn_daf, max_daf)
    daf_axis = range(1, max_daf)
    split_non_syn = defaultdict(list)
    if args.method == "SIFT":
        for sift_score, s, c in list_non_syn:
            assert (0 <= sift_score <= 1.0)
            if sift_score > 0.8:
                split_non_syn["pos"].append(c)
            elif sift_score > 0.3:
                split_non_syn["weak"].append(c)
            elif sift_score > 0.1:
                split_non_syn["neg"].append(c)
            else:
                split_non_syn["neg-strong"].append(c)
        dico_color = {"pos": RED, "weak": LIGHTGREEN, "neg": GREEN, "neg-strong": BLUE}
        dico_label = {"pos": "$0.8<SIFT$", "weak": "$0.3<SIFT<0.8$", "neg": "$0.1<SIFT<0.3$",
                      "neg-strong": "$SIFT<0.1$"}

        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        x = [sift_score for sift_score, s, c in list_non_syn if abs(s) < 20]
        y = [s for sift_score, s, c in list_non_syn if abs(s) < 20]
        plt.scatter(x, y, alpha=0.4, s=5.0)
        plt.xlabel("SIFT score")
        plt.ylabel("S given by Mutation-Selection")
        plt.tight_layout()
        plt.savefig(args.output.replace('.SIFT.pdf', '.SIFT_vs_MutSel.pdf'), format="pdf")
        plt.clf()
        plt.close("all")
    else:
        for s, c in list_non_syn:
            if s > 1.0:
                split_non_syn["pos"].append(c)
            elif s > 0.0:
                split_non_syn["pos-weak"].append(c)
            elif s > -1.0:
                split_non_syn["neg-weak"].append(c)
            elif s > -3.0:
                split_non_syn["neg"].append(c)
            else:
                split_non_syn["neg-strong"].append(c)
        dico_color = {"pos": RED, "pos-weak": YELLOW, "neg-weak": LIGHTGREEN, "neg": GREEN, "neg-strong": BLUE}
        dico_label = {"pos": "$1<S$", "pos-weak": "$0<S<1$", "neg-weak": "$-1<S<0$", "neg": "$-3<S<-1$",
                      "neg-strong": "$S<-3$"}

    for norm_syn in [True, False]:

        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        if norm_syn:
            plt.plot(daf_axis, [1.0] * len(daf_axis), linewidth=2.0, color="black")
        else:
            mean_count = np.mean(syn_count, axis=0)
            std_count = np.std(syn_count, axis=0)
            plt.scatter(daf_axis, mean_count, color="black")
            nb_snp = int(np.sum([1 for i in np.array(list_syn_daf).flatten() if i != 0]) / nbr_replicates)
            plt.plot(daf_axis, mean_count, linewidth=1.0, color="black",
                     label=f"Synonymous ({nb_snp} SNPs)")
            plt.fill_between(daf_axis, mean_count - std_count, mean_count + std_count, linewidth=1.0, color="black",
                             alpha=0.2)
        for key in dico_label:
            if len(split_non_syn[key]) == 0:
                continue
            non_syn_count = daf_to_count(split_non_syn[key], max_daf)

            if norm_syn:
                non_syn_count /= syn_count

            mean_count = np.mean(non_syn_count, axis=0)
            std_count = np.std(non_syn_count, axis=0)
            nb_snp = int(np.sum([1 for i in np.array(split_non_syn[key]).flatten() if i != 0]) / nbr_replicates)
            label = dico_label[key] + f" ({nb_snp} SNPs)"
            plt.scatter(daf_axis, mean_count, color=dico_color[key])
            plt.plot(daf_axis, mean_count, label=label, color=dico_color[key], linewidth=1.0)
            plt.fill_between(daf_axis, mean_count - std_count, mean_count + std_count, linewidth=1.0,
                             color=dico_color[key], alpha=0.2)

        plt.legend()
        plt.xlabel("Derived allele count")
        if norm_syn:
            plt.ylabel('Frequency of SNPs (relative to synonymous)')
        else:
            plt.ylabel("Frequency of SNPs")
        plt.tight_layout()
        output = args.output.replace('.pdf', '.normalize.pdf') if norm_syn else args.output
        plt.savefig(output, format="pdf")
        plt.clf()
        plt.close("all")
