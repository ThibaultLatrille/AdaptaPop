import argparse
import statsmodels.api as sm
import gzip
import numpy as np
from libraries import RED, GREEN, BLUE, LIGHTGREEN, YELLOW
from collections import defaultdict
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 256

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=False, type=str, dest="vcf", help="Input vcf file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output pdf file")
    args = parser.parse_args()

    subsample = 16
    vcf_file = gzip.open(args.vcf, 'rt')
    list_non_syn, list_syn_daf = [], []
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                header = {k: i for i, k in enumerate(vcf_line.strip().split("\t"))}
            continue

        info = vcf_line.strip().split("\t")[header["INFO"]]
        dico_info = {k: v for k, v in [s.split("=") for s in info.split(";") if "=" in s]}
        n = int(dico_info["SAMPLE_SIZE"])
        k = int(dico_info["COUNT"])
        count = np.random.hypergeometric(k, n - k, subsample) if n > subsample else k
        if count == 0 or count == min(n, subsample):
            continue
        if dico_info["SNP_TYPE"] == "Syn":
            list_syn_daf.append(count)
        elif dico_info["SNP_TYPE"] == "NonSyn":
            list_non_syn.append((float(dico_info["SEL_COEFF"]), count))

    vcf_file.close()

    syn_bincount = np.array(np.bincount(list_syn_daf)[1:]) / len(list_syn_daf)
    syn_counts = range(1, np.amax(list_syn_daf) + 1)
    split_non_syn = defaultdict(list)
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
            plt.plot(syn_counts, [1.0] * len(syn_counts), linewidth=2.0, color="black")
        else:
            plt.scatter(syn_counts, syn_bincount, color="black")
            plt.plot(syn_counts, syn_bincount, linewidth=1.0, color="black",
                     label=f"Synonymous ({len(list_syn_daf)} SNPs)")
        for key in dico_label:
            non_syn_bins_daf = split_non_syn[key]
            non_syn_bincount = np.array(np.bincount(non_syn_bins_daf)[1:]) / len(non_syn_bins_daf)
            if norm_syn:
                non_syn_bincount /= syn_bincount
            non_syn_counts = range(1, np.amax(non_syn_bins_daf) + 1)
            label = dico_label[key] + f" ({len(non_syn_bins_daf)} SNPs)"
            plt.scatter(non_syn_counts, non_syn_bincount, color=dico_color[key])
            plt.plot(non_syn_counts, non_syn_bincount, label=label, color=dico_color[key], linewidth=1.0)

        plt.legend()
        plt.xlabel("Derived allele count")
        if norm_syn:
            plt.ylabel('Frequency of SNPs (relative to synonymous)')
        else:
            plt.ylabel("Frequency of SNPs")
        plt.tight_layout()
        output = args.output.replace('.pdf', '.normalize.pdf') if norm_syn else args.output
        plt.savefig(output, format="pdf")
        plt.savefig(output.replace(".pdf", ".png"), format="png")
        plt.clf()
        plt.close("all")

    '''
    x = [s[0] for s in list_pairs]
    y = [s[1] for s in list_pairs]
    plt.scatter(x, y)
    model = sm.OLS(y, sm.add_constant(x))
    results = model.fit()
    b, a = results.params[0:2]
    idf = np.linspace(min(x), max(x), 100)
    plt.plot(idf, a * idf + b, '-',
             label=f"$y={a:.2g}x {'+' if float(b) > 0 else '-'} {abs(b):.2g}$ ($r^2={results.rsquared:.2g})$")
    '''
