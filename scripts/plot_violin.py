import argparse
import seaborn as sns
import pandas as pd
from collections import defaultdict
from libraries import format_pop, sp_to_color, RED, GREEN, sp_sorted
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

fontsize = 16
fontsize_legend = 14
my_dpi = 256
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output folder")
    args = parser.parse_args()

    pop2sp = {}
    groups = defaultdict(list)
    for filepath in args.tsv:
        _, model, sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")
        sp, pop, granu = filepath.split("/")[-2].replace("_", " ").split("-")
        pop2sp[format_pop(pop)] = sp
        groups[(model, sfs, granu)].append(filepath)

    order = [pop for pop, sp in sorted(pop2sp.items(), key=lambda kv: sp_sorted(*kv))]
    for (model, sfs, granu), files in groups.items():
        dico_nn, dico_ada = defaultdict(list), defaultdict(list)
        for filepath in files:
            df = pd.read_csv(filepath, sep="\t")
            ddf_nn = df[(~df["ADAPTIVE"]) & (df["OMEGA_A"] != "None")]
            ddf_ada = df[(df["ADAPTIVE"]) & (df["OMEGA_A"] != "None")]
            if len(ddf_nn) == 0 or len(ddf_ada) == 0:
                continue
            _, pop, _ = filepath.split("/")[-2].replace("_", " ").split("-")

            dico_ada["wA"].extend(ddf_ada["OMEGA_A"])
            dico_ada["pop"].extend([format_pop(pop)] * len(ddf_ada))
            dico_nn["wA"].extend(ddf_nn["OMEGA_A"])
            dico_nn["pop"].extend([format_pop(pop)] * len(ddf_nn))

        df_nn = pd.DataFrame(dico_nn)
        df_ada = pd.DataFrame(dico_ada)
        sub_pop = set(df_nn["pop"])
        nbr_pop = len(sub_pop)
        sub_order = [pop for pop in order if pop in sub_pop]

        fig, ax = plt.subplots(figsize=(210 * (nbr_pop + 1) / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        graph = sns.violinplot(x="pop", y="wA", data=df_nn, ax=ax, color=GREEN, inner="quartile", order=sub_order,
                               label=f"Nearly-neutral {granu}s (subsampled)")
        sns.stripplot(x="pop", y="wA", data=df_ada, ax=ax, color=RED, edgecolor="gray", size=10, jitter=0,
                      order=sub_order, label=f"Adaptive {granu}s")
        # plt.legend(fontsize=fontsize)

        ax.set_xlabel("", fontsize=fontsize_legend)
        ax.set_ylabel(r"$\omega_A$", fontsize=fontsize)
        # ax.set_title("{0} SFS at {1} level with {2}".format(sfs.capitalize(), granu, model), fontsize=fontsize)
        graph.axhline(0.0, color=GREEN, linestyle="--")
        if granu == "gene":
            graph.axhline(0.12, color=RED, linestyle="--")
        if granu == "site":
            graph.axhline(0.52, color=RED, linestyle="--")
        # plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
        for tick_label in graph.axes.get_xticklabels():
            tick_label.set_color(sp_to_color(pop2sp[tick_label.get_text()]))
            tick_label.set_fontsize(fontsize_legend)
        for tick_label in graph.axes.get_yticklabels():
            tick_label.set_fontsize(fontsize_legend)
        plt.tight_layout()
        plt.savefig(args.output + "/{0}-{1}-{2}.pdf".format(granu, sfs, model), format="pdf")
        plt.clf()
        plt.close("all")
