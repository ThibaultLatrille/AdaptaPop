import argparse
import seaborn as sns
import pandas as pd
from collections import defaultdict
from libraries import format_pop, sp_to_color, RED, GREEN
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt


def sp_gather(specie):
    # return specie + "-" if "Chlorocebus" in specie else "HS-OA-BT-"
    return ""


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
        groups[(model, sfs, granu, sp_gather(sp))].append(filepath)

    for (model, sfs, granu, sp), files in groups.items():
        dico_output = defaultdict(list)
        for filepath in files:
            df = pd.read_csv(filepath, sep="\t")
            if len(df["OMEGA_A"]) == 0 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) == 0:
                continue
            df = df[df["OMEGA_NA"] != "None"]
            _, pop, _ = filepath.split("/")[-2].replace("_", " ").split("-")

            dico_output["wA"].extend(df["OMEGA_A"])
            dico_output["regime"].extend(["adaptive" if i else "nearly-neutral" for i in df["ADAPTIVE"]])
            dico_output["pop"].extend([format_pop(pop)] * len(df))

        df = pd.DataFrame(dico_output)
        nbr_pop = len(set(df["pop"]))
        fig, ax = plt.subplots(figsize=(210 * (nbr_pop + 1) / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        palette = {"adaptive": "#EB6231", "nearly-neutral": "#8FB03E"}
        hue_order = ["nearly-neutral", "adaptive"]
        graph = sns.violinplot(x="pop", y="wA", hue="regime", data=df, ax=ax, palette=palette, hue_order=hue_order,
                               inner="quartile")
        plt.legend(fontsize=fontsize)
        '''
        graph = sns.violinplot(x="pop", y="wA", hue="regime", data=df, ax=ax, palette=palette,
                               hue_order=hue_order, inner=None, linewidth=0, saturation=0.4)
        sns.boxplot(x="pop", y="wA", hue="regime", data=df, palette=palette, boxprops={'zorder': 2}, ax=graph)
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles[2:], labels[2:], fontsize=fontsize)
        '''

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
        plt.savefig(args.output + "/{0}{1}-{2}-{3}.pdf".format(sp.replace(" ", "_"), granu, sfs, model), format="pdf")
        plt.clf()
        plt.close("all")
