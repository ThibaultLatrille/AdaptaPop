import argparse
import gzip
import numpy as np
import pandas as pd
from libraries_plot import plt, RED, GREEN, BLUE, LIGHTGREEN, YELLOW
from collections import defaultdict

my_dpi = 256
weak = ["A", "T"]
strong = ["G", "C"]
nbr_replicates = 10
subsample = 16


def daf_to_count(daf_list, min_n):
    array_daf = np.array(daf_list, dtype=np.int64).T
    return np.array([np.array(np.bincount(d, minlength=min_n)[1:]) / sum([1 for i in d if i != 0]) for d in array_daf])


omega_scaling = {"watterson": lambda i, n: 1.0 / i, "tajima": lambda i, n: n - i, "fay_wu": lambda i, n: i}


def theta(epsilon, daf_n, omega_scaling_method):
    transform = np.array([omega_scaling[omega_scaling_method](i, daf_n) for i in range(1, daf_n)])
    return sum(epsilon * transform) / sum(transform)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=False, type=str, dest="vcf", help="Input vcf file")
    parser.add_argument('-s', '--method', required=False, type=str, dest="method", help="Sel coeff parameter")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output pdf file")
    args = parser.parse_args()

    assert args.method in ["MutSel", "Omega", "Omega_0", "WS", "SW", "SIFT"]

    vcf_file = gzip.open(args.vcf, 'rt')
    list_non_syn = []
    max_daf_set = set()
    snp_cat = defaultdict(list)
    header = {}
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
        dico_info = {k: v for k, v in [s.split('=') for s in info.split(';') if '=' in s]}

        if args.method == "WS" and not ((dico_info["NUC_ANC"] in weak) and (dico_info["NUC_DER"] in strong)):
            continue

        if args.method == "SW" and not ((dico_info["NUC_ANC"] in strong) and (dico_info["NUC_DER"] in weak)):
            continue

        polarized = dico_info["POLARIZED"] == "True"
        if not polarized and args.method == "SIFT" and dico_info["SNP_TYPE"] == "NonSyn":
            continue

        sample_size = int(dico_info["SAMPLE_SIZE"])
        k = int(dico_info["COUNT_POLARIZED"])
        max_daf = min(sample_size, subsample)
        max_daf_set.add(max_daf)
        if sample_size > subsample:
            count = [np.random.hypergeometric(k, sample_size - k, subsample) for i in range(nbr_replicates)]
            count = [i if i != max_daf else 0 for i in count]
        else:
            count = [k]

        if np.all([i == 0 for i in count]):
            continue

        if dico_info["SNP_TYPE"] == "Syn":
            snp_cat["syn"].append(count)
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
                sift_find = sif_info.find('Sift=0')
                if sift_find == -1:
                    continue
                sift_score = float(sif_info[sift_find:].split("|")[2])
                list_non_syn.append((sift_score, S, count))
            elif np.isfinite(S):
                list_non_syn.append((S, count))
    vcf_file.close()

    if args.method == "SIFT":
        for sift_score, s, c in list_non_syn:
            assert (0 <= sift_score <= 1.0)
            if sift_score > 0.8:
                snp_cat["pos"].append(c)
            elif sift_score > 0.3:
                snp_cat["weak"].append(c)
            elif sift_score > 0.1:
                snp_cat["neg"].append(c)
            else:
                snp_cat["neg-strong"].append(c)
        dico_color = {"syn": "black", "pos": RED, "weak": LIGHTGREEN, "neg": GREEN, "neg-strong": BLUE}
        dico_label = {"syn": "Synonymous", "pos": "$0.8<SIFT$", "weak": "$0.3<SIFT<0.8$", "neg": "$0.1<SIFT<0.3$",
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
                snp_cat["pos"].append(c)
            elif s > 0.0:
                snp_cat["pos-weak"].append(c)
            elif s > -1.0:
                snp_cat["neg-weak"].append(c)
            elif s > -3.0:
                snp_cat["neg"].append(c)
            else:
                snp_cat["neg-strong"].append(c)
        dico_color = {"syn": "black", "pos": RED, "pos-weak": YELLOW, "neg-weak": LIGHTGREEN, "neg": GREEN,
                      "neg-strong": BLUE}
        dico_label = {"pos": "$1<S$", "pos-weak": "$0<S<1$", "syn": "Synonymous", "neg-weak": "$-1<S<0$",
                      "neg": "$-3<S<-1$", "neg-strong": "$S<-3$"}

    assert len(max_daf_set) == 1
    max_daf = max(max_daf_set)
    daf_axis = range(1, max_daf)
    scale = np.array([i for i in range(1, max_daf)])
    theta_dict = defaultdict(list)

    for scaled in [True, False]:
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for key in dico_label:
            if len(snp_cat[key]) == 0:
                continue

            non_syn_count = daf_to_count(snp_cat[key], max_daf)

            if scaled:
                non_syn_count *= scale

            mean_count = np.mean(non_syn_count, axis=0)
            std_count = np.std(non_syn_count, axis=0)
            nb_snp = int(np.sum([1 for i in np.array(snp_cat[key]).flatten() if i != 0]) / nbr_replicates)
            label = dico_label[key] + f" ({nb_snp} SNPs)"
            plt.scatter(daf_axis, mean_count, color=dico_color[key])
            plt.plot(daf_axis, mean_count, label=label, color=dico_color[key], linewidth=1.0)
            plt.fill_between(daf_axis, mean_count - std_count, mean_count + std_count, linewidth=1.0,
                             color=dico_color[key], alpha=0.2)
            if not scaled:
                theta_dict["category"].append(dico_label[key])
                for scaling_param in omega_scaling:
                    theta_dict[scaling_param].append(theta(mean_count, max_daf, scaling_param))

        plt.legend()
        plt.xlabel("Derived allele count")
        if scaled:
            plt.ylabel('Frequency of SNPs (scaled)')
        else:
            plt.ylabel("Frequency of SNPs")
        plt.tight_layout()
        output = args.output.replace('.pdf', '.normalize.pdf') if scaled else args.output
        plt.savefig(output, format="pdf")
        plt.clf()
        plt.close("all")

    df = pd.DataFrame(theta_dict)
    df.to_csv(args.output.replace('.pdf', '.tsv'), sep="\t", index=False)
