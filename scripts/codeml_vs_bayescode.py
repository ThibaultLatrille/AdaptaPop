import os
import argparse
import statsmodels.api as sm
from libraries import build_divergence_dico, filtered_table_omega, codontable
from libraries_plot import plt, np, LIGHTGREEN, BLUE

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


def open_fasta(fasta, open_all_sites=True):
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
        if ((column.count('X') + column.count('-')) == 0) or open_all_sites:
            sites_with_var.append(i)
            aa_list.append(seqs[0][i])
    return aa_list, sites_with_var


def vector_line(f, key, k):
    v_line = f.readline()
    assert v_line.startswith(f"{key}: ")
    v = [float(i) for i in v_line[3:].strip().replace("  ", " ").split(" ")]
    assert len(v) == k
    return v


def open_rst(file, aa_list):
    out_file = {}
    with open(file, "r") as f:
        line = f.readline()
        for line in f:
            if "MLEs of dN/dS (w) for site classes" in line:
                break
        if line == "":
            return out_file, 0

        k = int(line[line.rfind("(K=") + 3:line.rfind(")")])
        f.readline()
        p = vector_line(f, "p", k)
        w = vector_line(f, "w", k)
        gene_w = np.dot(p, w)

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
            if aa_list[i] != s[1]:
                s[-back] = np.nan
            i += 1
            out_file[int(s[0])] = float(s[-back])
    return out_file, gene_w


def all_sites(ctl):
    r = open(ctl, "r").read()
    return "cleandata = 0" in r


def build_codeml_dico(codeml_folder, list_ensg, template):
    print('Loading CODEML results.')
    site_codeml, gene_codeml = [], []
    filter_subset = {}
    for ensg in list_ensg:
        ctl = f"{codeml_folder}/{ensg}_NT_{template}/{ensg}_NT.ctl"
        fasta = f"{codeml_folder}/{ensg}_NT_{template}/{ensg}_NT.fasta"
        aa_list, sites_with_var = open_fasta(fasta, all_sites(ctl))
        rst = f"{codeml_folder}/{ensg}_NT_{template}/rst"
        if not os.path.exists(rst):
            continue
        out_file, gene_w = open_rst(rst, aa_list)
        if len(out_file) == 0:
            continue
        assert len(out_file) == len(sites_with_var)
        print(f"{ensg} with {template} model: {len(out_file)} sites")
        filter_subset[ensg] = sites_with_var
        site_codeml.extend(out_file.values())
        gene_codeml.append(gene_w)

    print(f'CODEML results loaded for {len(filter_subset)} genes.')
    return np.array(site_codeml), np.array(gene_codeml), filter_subset


def plot_codeml_bayescode(bayescode_w, codeml_w, output, level="sites"):
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    plt.subplot(1, 1, 1)
    plt.ylabel(f"BayesCode {level[:-1]} $\\omega$", fontsize=fontsize)
    plt.xlabel(f"CODEML {level[:-1]} $\\omega$", fontsize=fontsize)
    f = np.isfinite(codeml_w) & (bayescode_w[:, 1] < 10) & (bayescode_w[:, 1] > 0) & (codeml_w > 0)
    codeml_w = codeml_w[f]
    bayescode_w = bayescode_w[f, :]
    xmin, xmax = min(min(bayescode_w[:, 1]), min(codeml_w)), max(max(bayescode_w[:, 1]), max(codeml_w)) + 0.1
    plt.xlim((xmin, xmax))
    idf = np.linspace(xmin, xmax, 30)
    plt.plot(idf, idf, color="black")

    plt.errorbar(codeml_w, bayescode_w[:, 1],
                 yerr=[bayescode_w[:, 1] - bayescode_w[:, 0], bayescode_w[:, 2] - bayescode_w[:, 1]],
                 alpha=1.0, c=("#e6ebf3" if level == "sites" else "#adbfd9"), fmt='o', marker=None, mew=0,
                 zorder=-10, lw=(0.5 if level == "sites" else 1.0), markersize=0.0)
    plt.scatter(codeml_w, bayescode_w[:, 1], alpha=(0.4 if level == "sites" else 1.0), color=BLUE,
                label=f"${len(bayescode_w)}$ {level}", zorder=0, s=(1.5 if level == "sites" else 12.0))

    # outside = len([c for c, b in zip(codeml_w, bayescode_w) if c < b[0] or c > b[2]])
    # s = f"{outside} {level} out of {len(bayescode_w)} ({(100 * outside / len(bayescode_w)):.2f}%)"
    # s += f" are outside the {100 * (1 - float(pp) * 2):.0f}% confidence interval."
    # plt.title(s)
    if len(bayescode_w) > 1:
        model = sm.OLS(list(bayescode_w[:, 1]), codeml_w)
        results = model.fit()
        a = results.params[0]
        plt.plot(idf, a * idf, '-', label=f"$y={a:.3g}x$ ($r^2={results.rsquared:.3g})$", color=LIGHTGREEN)
    plt.legend(fontsize=fontsize_legend, loc="lower right")
    plt.ylim((xmin, xmax))
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(f"{output}.png", format="png")
    plt.clf()
    plt.close('all')


def main(codeml_folder, bayescode_folder, pp, output):
    list_codeml = sorted(os.listdir(codeml_folder))
    templates = sorted(set([i.split("_")[-1] for i in list_codeml]))
    for template in templates:
        list_ensg = [i.replace(f"_NT_{template}", "") for i in list_codeml if i.split("_")[-1] == template]
        site_codeml, gene_codeml, filter_subset = build_codeml_dico(codeml_folder, list_ensg, template)

        _, gene_omega_dico = build_divergence_dico(bayescode_folder, list_ensg, gene_level=True, pp=pp)
        gene_bayescode = filtered_table_omega(gene_omega_dico, filter_subset, gene_level=True)
        plot_codeml_bayescode(gene_bayescode, gene_codeml, f"{output}.gene.{template}", level="genes")

        _, site_omega_dico = build_divergence_dico(bayescode_folder, list_ensg, gene_level=False, pp=pp)
        site_bayescode = filtered_table_omega(site_omega_dico, filter_subset, gene_level=False)
        plot_codeml_bayescode(site_bayescode, site_codeml, f"{output}.site.{template}", level="sites")
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
