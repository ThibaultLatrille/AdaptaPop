import argparse
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import itertools
from libraries import tex_f, format_pop, sp_to_color, sp_sorted
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
    for tick_label in ax.axes.get_xticklabels():
        tick_label.set_color(sp_to_color(tick_label.get_text()))
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


def extend_pop(p, sample):
    fp = format_pop(p)
    if sample[p] == fp:
        return p
    else:
        return f"{sample[p]} ({fp})"


header = ["Species", "Population", "SFS", "Level", "Model",
          "$\\omega_{\\textrm{A}}^{\\textrm{S}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{S}}$",
          "$\\omega^{\\textrm{S}}$", "$\\alpha^{\\textrm{S}}$",
          "$\\omega_{\\textrm{A}}^{\\textrm{N}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{N}}$",
          "$\\omega^{\\textrm{N}}$", "$\\alpha^{\\textrm{N}}$",
          "p-value", "$\\pi_{\\textrm{N}}$"]

header = ["Species", "Population", "SFS", "Level", "Model",
          "$\\color{RED}{\\omega_{\\mathrm{A}}^{\\mathrm{AR}}}\\color{black}$",
          "$\\color{GREEN}{\\omega_{\\mathrm{A}}^{\\mathrm{NNR}}}\\color{black}$",
          "$\\Delta_{\\omega_{\\mathrm{A}}}$",
          "p-value", "$\\pi_{\\textrm{S}}$"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    parser.add_argument('-s', '--sample_list', required=False, type=str, dest="sample_list", help="Sample list file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df = df.groupby(["species", "pop", "sfs", "granularity", "model"]).max().reset_index()

    dico_matrix, dico_delta_wa = defaultdict(dict), defaultdict(dict)
    pop2sp = {}
    for sfs, granularity, model in itertools.product(["folded", "unfolded"], ["gene", "site"],
                                                     ["grapes", "polyDFE", "dfem", "MK", "aMK"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)]
        m = "{0}s - {1} SFS - {2}".format(granularity.capitalize(), sfs, model)
        for pop in ddf["pop"].values:
            pop_ddf = ddf[ddf["pop"] == pop]
            if len(pop_ddf["pop"]) != 0:
                assert len(pop_ddf["pop"]) == 1
                fpop = format_pop(pop)
                pop2sp[fpop] = pop_ddf["species"].values[0]
                dico_matrix[fpop][m] = pop_ddf["p_val"].values[0]
                dico_delta_wa[fpop][m] = pop_ddf["wA_Selected"].values[0] - pop_ddf["wA_Neutral"].values[0]

    models = set()
    for s in dico_matrix.values():
        models = models.union(s.keys())
    models = list(sorted(models))

    species = [pop for pop, sp in sorted(pop2sp.items(), key=lambda kv: sp_sorted(*kv))]
    p_matrix, delta_wa_matrix = np.ones((len(species), len(models))), np.ones((len(species), len(models)))
    p_matrix[:], delta_wa_matrix[:] = np.NaN, np.NaN

    for id_species, sp in enumerate(species):
        for id_model_set, model_set in enumerate(models):
            if sp in dico_matrix and model_set in dico_matrix[sp]:
                p_matrix[id_species, id_model_set] = dico_matrix[sp][model_set]
                delta_wa_matrix[id_species, id_model_set] = dico_delta_wa[sp][model_set]

    fig, ax = plt.subplots()

    YlGn = matplotlib.cm.YlGn
    im, cbar = heatmap(p_matrix.T, models, species, ax=ax, cmap=YlGn, cbarlabel=r"$p_{\mathrm{value}}$")
    texts = annotate_heatmap(im, valfmt=lambda p: "0" if abs(p) < 1e-1 else "{0:.1g}".format(p), fontsize=6)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".pval.png"), format="png")
    plt.savefig(args.output.replace(".tex", ".pval.pdf"), format="pdf")
    plt.clf()

    fig, ax = plt.subplots()
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

    df = df.iloc[df.apply(lambda r: sp_sorted(format_pop(r["pop"]), r["species"]), axis=1).argsort()]
    df["species"] = df.apply(lambda r: "Ovis orientalis" if r["pop"] == "IROO" else r["species"], axis=1)
    df["species"] = df.apply(lambda r: "Ovis vignei" if r["pop"] == "IROV" else r["species"], axis=1)
    df["species"] = df.apply(lambda r: "Capra aegagrus" if r["pop"] == "IRCA" else r["species"], axis=1)

    dico_sample = {r["SampleName"].replace("_", " "): r["Location"] for _, r in
                   list(pd.read_csv(args.sample_list, sep='\t').iterrows())}

    df["pop"] = df.apply(lambda r: extend_pop(r["pop"], dico_sample), axis=1)

    o = open(args.output, 'w')
    sub_header = [i for i in header if i not in ["SFS", "Level", "Model"]]
    for sfs, granularity, model in itertools.product(["folded", "unfolded"], ["gene", "site"],
                                                     ["grapes", "polyDFE", "dfem", "MK", "aMK"]):
        ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)].drop(
            ["sfs", "granularity", "model"], axis=1)
        if len(ddf["species"]) == 0: continue
        o.write("\\subsection{" + "{0} SFS at {1} level - {2}".format(sfs.capitalize(), granularity, model) + "} \n")
        o.write("\\begin{center}\n")
        o.write("\\includegraphics[width=\\linewidth]{ViolinPlot/" +
                "{0}-{1}-{2}.pdf".format(granularity, sfs, model) + "} \n")
        o.write("\\begin{adjustbox}{width = 1\\textwidth}\n")
        o.write(ddf.to_latex(index=False, escape=False, float_format=tex_f, header=sub_header))
        o.write("\\end{adjustbox}\n")
        o.write("\\end{center}\n")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {0}/main-table.tex".format(
        os.path.dirname(args.output))
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
