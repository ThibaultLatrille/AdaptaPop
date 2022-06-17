import argparse
import os
from collections import defaultdict
from functools import reduce
from libraries_plot import *
from matplotlib import cm

header = {"pop": "Population", "species": "Species",
          "wA_Test": "\\specialcell{$\\omega_{\\mathrm{A}}$ \\\\ Adaptive}",
          "wA_Control": "\\specialcell{$\\left< \\omega_{\\mathrm{A}} \\right>$ \\\\ Nearly-neutral}",
          "wA_Delta": "$\\Delta \\omega_{\\mathrm{A}} $",
          "w_Test": "\\specialcell{$d_{\\mathrm{N}} / d_{\\mathrm{S}}$ \\\\ Adaptive}",
          "w_Control": "\\specialcell{$\\left< d_{\\mathrm{N}} / d_{\\mathrm{S}} \\right>$ \\\\ Nearly-neutral}",
          "r": "$\\frac{\\Delta\\omega_{\\mathrm{A}}}{\\Delta\\omega_{\\mathrm{A}}^{\\mathrm{phy}}}$",
          "wA_pval": "$p_{\\mathrm{v}}$", "w_pval": "$p_{\\mathrm{v}}$",
          "wA_pval_adj": "$p_{\\mathrm{v}}^{\\mathrm{adj}}$", "w_pval_adj": "$p_{\\mathrm{v}}^{\\mathrm{adj}}$",
          "P_S": "$\\pi_{\\textrm{S}}$"}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    parser.add_argument('-b', '--bounds', required=False, type=str, dest="bounds", default="div",
                        help="Integral inferior bound for the calculation of alpha by polyDFE")
    parser.add_argument('-s', '--sample_list', required=False, type=str, dest="sample_list", help="Sample list file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df["wA_Delta"] = df["wA_Test"] - df["wA_Control"]
    df["w_Delta"] = df["Δw_Test"] - df["Δw_Control"]
    df["r"] = df["wA_Delta"] / df["w_Delta"]

    dico_pval, dico_delta_wa = defaultdict(dict), defaultdict(dict)
    pop2sp = {}
    for (sfs, level, method, pp, model), ddf in df.groupby(["sfs", "level", "method", "pp", "model"]):
        m = f"{method} {level}s (pp={pp}) - {sfs[0]}SFS - {model}"
        for pop, pop_ddf in ddf.groupby(["pop"]):
            assert len(pop_ddf["pop"]) == 1
            fpop = format_pop(pop)
            pop2sp[fpop] = pop_ddf["species"].values[0]
            dico_pval[fpop][m] = pop_ddf["wA_pval"].values[0]
            dico_delta_wa[fpop][m] = pop_ddf["wA_Delta"].values[0]

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

    YlGn = cm.YlGn
    im, cbar = heatmap(pval_matrix.T, models, species, ax=ax, cmap=YlGn, cbarlabel=header["w_pval"])
    texts = annotate_heatmap(im, valfmt=lambda p: "0" if abs(p) < 1e-1 else "{0:.1g}".format(p), fontsize=6)
    plt.tight_layout()
    plt.savefig(args.output.replace(".tex", ".pval.pdf"), format="pdf")
    plt.clf()

    fig, ax = plt.subplots(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    RdBu = cm.get_cmap('RdBu_r')
    start = np.nanmin(delta_wa_matrix)
    midpoint = - start / (np.nanmax(delta_wa_matrix) - start)
    shifted_RdBu = shiftedColorMap(RdBu, midpoint=midpoint, name='shifted')
    im, cbar = heatmap(delta_wa_matrix.T, models, species, ax=ax, cmap=shifted_RdBu,
                       cbarlabel=header["wA_Delta"])
    texts = annotate_heatmap(im, valfmt=lambda p: "{0:.2f}".format(p), div=True, fontsize=5)
    plt.tight_layout()
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

    for (pp, sfs, level, method, model), ddf in df.groupby(["pp", "sfs", "level", "method", "model"]):
        assert len(ddf["species"]) != 0

        o.write(
            "\\subsection{" + f"{method} at {level} level (pp={pp}) - {sfs[0]}SFS - {model.replace('_', ' ')}" + "} \n")
        o.write("\\begin{center}\n")

        for prefix in ["wA", "w"]:
            if prefix == "wA":
                columns = ["pop", "species", "wA_Test", "wA_Control", "wA_Delta", "wA_pval", "wA_pval_adj", "r", "P_S"]
            else:
                columns = ["pop", "species", "w_Test", "w_Control", "w_pval", "w_pval_adj", "P_S"]
            ddf = adjusted_holm_pval(ddf, prefix=prefix + "_")
            o.write("\\includegraphics[width=\\linewidth]{ViolinPlot-" +
                    f"{args.bounds}/{level}-{method}-{pp}-{sfs}-{model}-{prefix}.pdf" + "} \n")
            o.write("\\begin{adjustbox}{width = 1\\textwidth}\n")
            o.write(ddf.to_latex(index=False, escape=False, float_format=tex_f, column_format=column_format(columns),
                                 header=[header[i] for i in columns], columns=columns))
            o.write("\\end{adjustbox}\n")
            o.write("\\newpage\n")
        o.write("\\end{center}\n")

    for (sfs, model, pp), ddf in df.groupby(["sfs", "model", "pp"]):
        ddg = {}
        for (level, method), ddf_g in ddf.groupby(["level", "method"]):
            adjusted_holm_pval(ddf_g, prefix="wA_")
            ddg[f"{level}-{method}"] = ddf_g

        list_df = list(ddg.values())
        merge = reduce(lambda l, r: pd.merge(l, r, how="inner", on=["pop", "species"]), list_df)

        columns, df_columns = ["pop", "species"], ["wA_Delta", "wA_pval_adj", "r"]
        sub_header = [header[i] for i in columns]

        columns.append("P_S_x")
        sub_header.append(header["P_S"])

        for key in ["_x", "_y", ""][:len(list_df)]:
            columns += [f"{i}{key}" for i in df_columns]
            sub_header += [header[i] for i in df_columns]
        o.write("\\subsection{" + f"{sfs[0]}SFS - {model.replace('_', ' ')} (pp={pp})" + "} \n")
        o.write("\\begin{center}\n")
        o.write("\\begin{adjustbox}{width = 1\\textwidth}\n")
        o.write(merge.to_latex(index=False, escape=False, float_format=tex_f, column_format=column_format(merge),
                               header=sub_header, columns=columns))
        o.write("\\end{adjustbox}\n")
        o.write("\\end{center}\n")
        o.write("\\newpage\n")
    o.close()

    folder = os.path.dirname(args.output)
    tex = f"{os.path.dirname(args.output)}/main-table-{args.bounds}.tex"
    os.system(f"sed 's/results.tex/results-{args.bounds}.tex/g' {folder}/main-table.tex > {tex}")
    tex_to_pdf = f"pdflatex -synctex=1 -interaction=nonstopmode -output-directory={folder} {tex}"
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
