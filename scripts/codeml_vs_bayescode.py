import os
import argparse
import statsmodels.api as sm
from libraries import build_divergence_dico, filtered_table_omega, codontable
from libraries_plot import plt, np, RED_RGB

error_kwargs = {"lw": .5, "zorder": -1}
my_dpi = 256
fontsize = 14
fontsize_legend = 12


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return f"{s}"


def translate_aa(seq):
    return [codontable[seq[i:i + 3]] for i in range(0, len(seq), 3)]


def open_fasta(fasta, all_sites=True):
    # Open a fasta file and return a list of sites with variation
    f = open(fasta, "r")
    seqs, aa_list, sites_with_var = [], [], []
    for line in f:
        if line[0] == ">":
            continue
        seqs.append(translate_aa(line.strip()))
    f.close()
    assert len(set([len(s) for s in seqs])) == 1
    for i in range(len(seqs[0])):
        column = "".join([s[i] for s in seqs])
        if ((column.count('X') + column.count('-')) == 0) or all_sites:
            sites_with_var.append(i)
            aa_list.append(seqs[0][i])
    return aa_list, sites_with_var


def open_rst(file, aa_list):
    out_file = {}
    with open(file, "r") as f:
        line = f.readline()
        for line in f:
            if "Naive Empirical Bayes (NEB) probabilities " in line:
                break

        back = line[line.rfind("postmean_w"):].count("&") + 1
        f.readline()
        f.readline()
        i = 0
        for line in f:
            s = line.replace("  ", " ").strip().split(" ")
            if len(s) == 1:
                break
            assert s[1] == aa_list[i]
            i += 1
            out_file[int(s[0])] = float(s[-back])

    return out_file


def all_sites(ctl):
    r = open(ctl, "r").read()
    return "cleandata = 0" in r


def build_codeml_dico(codeml_folder, list_ensg, template):
    print('Loading CODEML results.')
    codeml_w = []
    filter_subset = {}
    for ensg in list_ensg:
        ctl = f"{codeml_folder}/{ensg}_NT_{template}/{ensg}_NT.ctl"
        fasta = f"{codeml_folder}/{ensg}_NT_{template}/{ensg}_NT.fasta"
        aa_list, sites_with_var = open_fasta(fasta, all_sites(ctl))
        rst = f"{codeml_folder}/{ensg}_NT_{template}/rst"
        out_file = open_rst(rst, aa_list)
        if len(out_file) == 0:
            continue
        assert len(out_file) == len(sites_with_var)
        print(f"{ensg} with model {template}: {len(out_file)} sites")
        filter_subset[ensg] = sites_with_var
        codeml_w.extend(out_file.values())

    print('CODEML results loaded.')
    return codeml_w, filter_subset


def main(codeml_folder, bayescode_folder, pp, output):
    list_codeml = sorted(os.listdir(codeml_folder))
    templates = sorted(set([i.split("_")[-1] for i in list_codeml]))
    for template in templates:
        list_ensg = [i.replace(f"_NT_{template}", "") for i in list_codeml if i.split("_")[-1] == template]
        codeml_w, filter_subset = build_codeml_dico(codeml_folder, list_ensg, template)
        _, dico_omega = build_divergence_dico(bayescode_folder, list_ensg, gene_level=False, pp=pp)
        bayescode_w = filtered_table_omega(dico_omega, filter_subset, gene_level=False)

        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        plt.subplot(1, 1, 1)
        plt.ylabel("BayesCode $\\omega$", fontsize=fontsize)
        plt.xlabel("CODEML $\\omega$", fontsize=fontsize)
        xmin, xmax = 0, max(max(bayescode_w[:, 1]), max(codeml_w)) + 0.1
        plt.xlim((xmin, xmax))
        idf = np.linspace(xmin, xmax, 30)
        plt.plot(idf, idf, color="black")
        plt.errorbar(codeml_w, bayescode_w[:, 1],
                     yerr=[bayescode_w[:, 1] - bayescode_w[:, 0], bayescode_w[:, 2] - bayescode_w[:, 1]],
                     alpha=0.4, label=r"${0}$ sites".format(len(bayescode_w)), c=RED_RGB,
                     fmt='o', marker=None, mew=0, ecolor=RED_RGB, zorder=10, lw=.5, markersize=3.0)
        # Regression line
        model = sm.OLS(list(bayescode_w[:, 1]), sm.add_constant(codeml_w))
        results = model.fit()
        b, a = results.params[0:2]
        plt.plot(idf, a * idf + b, 'r-',
                 label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(float(a), abs(float(b)), results.rsquared,
                                                                          "+" if float(b) > 0 else "-"))
        plt.legend(fontsize=fontsize_legend, loc="lower right")
        plt.ylim((xmin, xmax))
        plt.xticks(fontsize=fontsize_legend)
        plt.tight_layout()
        plt.savefig(f"{output}.{template}.pdf", format="pdf")
        plt.clf()
        plt.close('all')

    print('Plot completed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bayescode_folder', required=True, type=str, dest="bayescode",
                        help="folder containing BayesCode results")
    parser.add_argument('-c', '--codeml_folder', required=True, type=str, dest="codeml",
                        help="folder containing CODEML results")
    parser.add_argument('-p', '--pp', required=True, type=str, dest="pp", help="Posterior probability")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.codeml, args.bayescode, args.pp, args.output)
