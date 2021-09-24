import argparse
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import itertools
from libraries import tex_f
import matplotlib

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    im = ax.imshow(data, **kwargs)
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw, orientation='horizontal')
    cbar.ax.set_xlabel(cbarlabel)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
             rotation_mode="anchor")
    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
    threshold = im.norm(data.max()) / 2.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j]), **kw)
            texts.append(text)
    return texts


def format_pop(t):
    if " " in t:
        return "".join([s[0] for s in t.split(" ")])
    else:
        return t


def format_pval(p):
    if abs(p) < 1e-1:
        return "0"
    else:
        return "{0:.1g}".format(p)


header = ["Species", "Population", "SFS", "Level", "Model",
          "$\\omega_{\\textrm{A}}^{\\textrm{S}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{S}}$",
          "$\\omega^{\\textrm{S}}$", "$\\alpha^{\\textrm{S}}$",
          "$\\omega_{\\textrm{A}}^{\\textrm{N}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{N}}$",
          "$\\omega^{\\textrm{N}}$", "$\\alpha^{\\textrm{N}}$",
          "p-value"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df = df.groupby(["species", "pop", "sfs", "granularity", "model"]).max().reset_index()
    df.sort_values(by=["species", "pop", "granularity", "sfs", "model"], inplace=True)

    dico_matrix, dico_delta_wa = defaultdict(dict), defaultdict(dict)
    for sfs, granularity, model in itertools.product(["folded", "unfolded"], ["gene", "site"],
                                                     ["grapes", "polyDFE", "dfem"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)]
        m = "{0}s with {1} on {2} SFS".format(granularity.capitalize(), model, sfs)
        for pop in ddf["pop"].values:
            pop_ddf = ddf[ddf["pop"] == pop]
            if len(pop_ddf["pop"]) != 0:
                assert len(pop_ddf["pop"]) == 1
                sp = "{0}-{1}".format(pop_ddf["species"].values[0].split(" ")[0],
                                      format_pop(pop_ddf["pop"].values[0]))
                dico_matrix[sp][m] = pop_ddf["p_val"].values[0]
                dico_delta_wa[sp][m] = pop_ddf["wA_Selected"].values[0] - pop_ddf["wA_Neutral"].values[0]

    models = set()
    for s in dico_matrix.values():
        models = models.union(s.keys())
    models = list(sorted(models))
    species = [str(i) for i in dico_matrix.keys()]
    p_matrix, delta_wa_matrix = np.ones((len(species), len(models))), np.ones((len(species), len(models)))
    p_matrix[:], delta_wa_matrix[:] = np.NaN, np.NaN

    for id_species, sp in enumerate(species):
        for id_model_set, model_set in enumerate(models):
            p_matrix[id_species, id_model_set] = dico_matrix[sp][model_set]
            delta_wa_matrix[id_species, id_model_set] = dico_delta_wa[sp][model_set]

    fig, ax = plt.subplots()
    im, cbar = heatmap(p_matrix.T, models, species, ax=ax, cmap="YlGn", cbarlabel=r"$p_{\mathrm{value}}$")
    texts = annotate_heatmap(im, valfmt=format_pval, fontsize=6)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".pval.png"), format="png")
    plt.savefig(args.output.replace(".tex", ".pval.pdf"), format="pdf")
    plt.clf()

    fig, ax = plt.subplots()
    im, cbar = heatmap(delta_wa_matrix.T, models, species, ax=ax, cmap="YlGn",
                       cbarlabel=r"$\Delta_{\omega_{\mathrm{A}}}$")
    texts = annotate_heatmap(im, valfmt=format_pval, fontsize=6)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".delta_wa.png"), format="png")
    plt.savefig(args.output.replace(".tex", ".delta_wa.pdf"), format="pdf")
    plt.clf()
    plt.close('all')

    o = open(args.output, 'w')
    # o.write(df.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=header))
    # o.write("\\newpage\n")

    # for model in ["grapes", "polyDFE", "dfem"]:
    #     ddf = df[df["model"] == model].drop(["model"], axis=1)
    #     if len(ddf["species"]) == 0: continue
    #     o.write("\\subsection{" + model + "} \n")
    #     o.write(ddf.to_latex(index=False, escape=False, longtable=True, float_format=tex_f,
    #                          header=[i for i in header if i != "Model"]))
    #     o.write("\\newpage\n")

    sub_header = [i for i in header if i not in ["SFS", "Level", "Model"]]
    for sfs, granularity, model in itertools.product(["folded", "unfolded"], ["gene", "site"],
                                                     ["grapes", "polyDFE", "dfem"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)].drop(
            ["sfs", "granularity", "model"], axis=1)
        if len(ddf["species"]) == 0: continue
        o.write("\\subsection{" + "{0} SFS at {1} level - {2}".format(sfs.capitalize(), granularity, model) + "} \n")
        o.write(ddf.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=sub_header))
        # o.write("\\newpage\n")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {0}/main-table.tex".format(
        os.path.dirname(args.output))
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
