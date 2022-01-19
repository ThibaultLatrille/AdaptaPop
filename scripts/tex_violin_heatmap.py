import argparse
import pandas as pd
import os
import numpy as np
from collections import defaultdict
from itertools import product
from libraries import tex_f, format_pop, sp_sorted, sort_df, adjusted_holm_pval
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
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


def annotate_heatmap(im, data=None, div=False, valfmt="{x:.2f}", textcolors=("white", "black"), **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
    threshold_high = im.norm(data.max()) * (0.75 if div else 0.5)
    threshold_low = (im.norm(data.max()) * 0.25) if div else 0.0
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(threshold_high >= im.norm(data[i, j]) >= threshold_low)])
            text = im.axes.text(j, i, valfmt(data[i, j]), **kw)
            texts.append(text)
    return texts


def format_pval(p):
    return "0" if abs(p) < 1e-1 else "{0:.1g}".format(p)


def format_domega(p):
    return "{0:.2g}".format(p)


def column_format(data_frame):
    return "|" + "|".join(["l"] * 2 + ["r"] * (len(data_frame) - 2)) + "|"

my_dpi = 256
omega_a = "\\specialcell{$\\omega_{\\mathrm{A}}$ \\\\ test}"
omega_a_mean = "\\specialcell{$\\left< \\omega_{\\mathrm{A}} \\right>$ \\\\ control}"
omega = "\\specialcell{$d_{\\mathrm{N}} / d_{\\mathrm{S}}$ \\\\ test}"
omega_mean = "\\specialcell{$\\left< d_{\\mathrm{N}} / d_{\\mathrm{S}} \\right>$ \\\\ control}"
chi = "$\\chi$"
pv = "$p_{\\mathrm{v}}$"
pv_adj = "$p_{\\mathrm{v}}^{\\mathrm{adj}}$"
pi_s = "$\\pi_{\\textrm{S}}$"

header_wa = ["Species", "Population", "SFS", "Level", "Model", omega_a, omega_a_mean, pv, pv_adj, chi, pi_s]
header_w = ["Species", "Population", "SFS", "Level", "Model", omega, omega_mean, pv, pv_adj, pi_s]
header_main = ["Species", "Population", omega_a, omega_a_mean, pv, pv_adj, chi, omega_a, omega_a_mean, pv, pv_adj, chi, pi_s]
invert_prefix = {"w": "wA", "wA": "w"}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    parser.add_argument('-s', '--sample_list', required=False, type=str, dest="sample_list", help="Sample list file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df["chi"] = (df["wA_Test"] - df["wA_Control"]) / (df["Δw_Test"] - df["Δw_Control"])
    df.pop("Δw_Test")
    df.pop("Δw_Control")
    df = df.groupby(["species", "pop", "sfs", "granularity", "model"]).max().reset_index()

    dico_pval, dico_delta_wa = defaultdict(dict), defaultdict(dict)
    pop2sp = {}
    for sfs, granularity, model in product(["folded", "unfolded"], ["gene", "site", "site strong"],
                                           ["grapes", "polyDFE", "dfem", "MK", "aMK"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)]
        m = "{0}s - {1}SFS - {2}".format(granularity.capitalize(), sfs[0], model)
        for pop in ddf["pop"].values:
            pop_ddf = ddf[ddf["pop"] == pop]
            if len(pop_ddf["pop"]) != 0:
                assert len(pop_ddf["pop"]) == 1
                fpop = format_pop(pop)
                pop2sp[fpop] = pop_ddf["species"].values[0]
                dico_pval[fpop][m] = pop_ddf["wA_pval"].values[0]
                dico_delta_wa[fpop][m] = pop_ddf["wA_Test"].values[0] - pop_ddf["wA_Control"].values[0]

    models = set()
    for s in dico_pval.values():
        models = models.union(s.keys())
    models = list(sorted(models))

    species = [pop for pop, sp in sorted(pop2sp.items(), key=lambda kv: sp_sorted(*kv))]
    pval_matrix, delta_wa_matrix = np.ones((len(species), len(models))), np.ones((len(species), len(models)))
    pval_matrix[:], delta_wa_matrix[:] = np.NaN, np.NaN

    for id_species, sp in enumerate(species):
        for id_model_set, model_set in enumerate(models):
            if sp in dico_pval and model_set in dico_pval[sp]:
                pval_matrix[id_species, id_model_set] = dico_pval[sp][model_set]
                delta_wa_matrix[id_species, id_model_set] = dico_delta_wa[sp][model_set]

    fig, ax = plt.subplots(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)

    YlGn = matplotlib.cm.YlGn
    im, cbar = heatmap(pval_matrix.T, models, species, ax=ax, cmap=YlGn, cbarlabel=r"$p_{\mathrm{value}}$")
    texts = annotate_heatmap(im, valfmt=lambda p: "0" if abs(p) < 1e-1 else "{0:.1g}".format(p), fontsize=6)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".pval.png"), format="png")
    plt.savefig(args.output.replace(".tex", ".pval.pdf"), format="pdf")
    plt.clf()

    fig, ax = plt.subplots(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    RdBu = matplotlib.cm.get_cmap('RdBu_r')
    start = np.nanmin(delta_wa_matrix)
    midpoint = - start / (np.nanmax(delta_wa_matrix) - start)
    shifted_RdBu = shiftedColorMap(RdBu, midpoint=midpoint, name='shifted')
    im, cbar = heatmap(delta_wa_matrix.T, models, species, ax=ax, cmap=shifted_RdBu,
                       cbarlabel="$\\Delta_{\\omega_{\\mathrm{A}}}$")
    texts = annotate_heatmap(im, valfmt=lambda p: "{0:.2f}".format(p), div=True, fontsize=5)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".delta_wa.png"), format="png")
    plt.savefig(args.output.replace(".tex", ".delta_wa.pdf"), format="pdf")
    plt.clf()
    plt.close('all')

    df = sort_df(df, args.sample_list)

    o = open(args.output, 'w')
    o.write("\\subsection{Heatmap} \n")
    o.write("\\begin{center}\n")
    o.write("\\includegraphics[width=\\linewidth]{" + args.output.replace(".tex", ".pval.pdf") + "} \n")
    o.write("\\includegraphics[width=\\linewidth]{" + args.output.replace(".tex", ".delta_wa.pdf") + "} \n")
    o.write("\\end{center}\n")

    for sfs, granularity, model in product(["folded", "unfolded"], ["gene", "site", "site strong"],
                                           ["grapes", "polyDFE", "dfem", "MK", "aMK"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)].drop(
            ["sfs", "granularity", "model"], axis=1)
        if len(ddf["species"]) == 0:
            continue

        o.write("\\subsection{" + "{0} SFS at {1} level - {2}".format(sfs.capitalize(), granularity, model) + "} \n")
        o.write("\\begin{center}\n")

        for prefix in ["wA", "w"]:
            sub_header = [i for i in (header_wa if prefix == "wA" else header_w) if i not in ["SFS", "Level", "Model"]]
            subdf = ddf.copy()
            subdf = adjusted_holm_pval(subdf, prefix=prefix + "_")
            subdf.pop(f"{invert_prefix[prefix]}_pval")
            subdf.pop(f"{invert_prefix[prefix]}_Test")
            subdf.pop(f"{invert_prefix[prefix]}_Control")
            if prefix == "w":
                subdf.pop("chi")
            else:
                subdf["chi"] = subdf.pop("chi")
            subdf["P_S"] = subdf.pop("P_S")
            o.write("\\includegraphics[width=\\linewidth]{ViolinPlot/" +
                    "{0}-{1}-{2}-{3}.pdf".format(granularity, sfs, model, prefix) + "} \n")
            o.write("\\begin{adjustbox}{width = 1\\textwidth}\n")
            o.write(subdf.to_latex(index=False, escape=False, float_format=tex_f, column_format=column_format(subdf),
                                   header=sub_header))
            o.write("\\end{adjustbox}\n")
            o.write("\\newpage\n")
        o.write("\\end{center}\n")

    sub_header_main = [i for i in header_main if i not in ["SFS", "Level", "Model"]]
    for sfs, model in product(["folded", "unfolded"], ["grapes", "polyDFE", "dfem", "MK", "aMK"]):
        ddf = df[(df["sfs"] == sfs) & (df["model"] == model)].drop(["sfs", "model"], axis=1)
        if len(ddf["species"]) == 0:
            continue

        for pop_col in ["w_pval", "w_Test", "w_Control"]:
            ddf.pop(pop_col)

        ddg = []
        for granularity in ["gene", "site"]:
            ddf_g = ddf[ddf["granularity"] == granularity].drop(["granularity"], axis=1)
            # ps = ddf_g.pop("P_S")
            adjusted_holm_pval(ddf_g, prefix="wA_")
            ddf_g["chi"] = ddf_g.pop("chi")
            if granularity == "gene":
                ddf_g.pop("P_S")
            else:
                ddf_g["P_S"] = ddf_g.pop("P_S")
            ddg.append(ddf_g)

        merge = pd.merge(ddg[0], ddg[1], how="inner", on=["pop", "species"])
        o.write("\\subsection{" + "{0} SFS - {1}".format(sfs.capitalize(), model) + "} \n")
        o.write("\\begin{center}\n")
        o.write("\\begin{adjustbox}{width = 1\\textwidth}\n")
        o.write(merge.to_latex(index=False, escape=False, float_format=tex_f, column_format=column_format(merge),
                               header=sub_header_main))
        o.write("\\end{adjustbox}\n")
        o.write("\\end{center}\n")
        o.write("\\newpage\n")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {0}/main-table.tex".format(
        os.path.dirname(args.output))
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
