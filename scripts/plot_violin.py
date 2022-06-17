import argparse
import seaborn as sns
from collections import defaultdict
from libraries_plot import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output folder")
    args = parser.parse_args()

    pop2sp = {}
    groups = defaultdict(list)
    for filepath in args.tsv:
        _, model, sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")[:3]
        sp, pop, level, method, pp = filepath.split("/")[-2].replace("_", " ").split("-")[:5]
        pop2sp[format_pop(pop)] = sp
        groups[(model, sfs, level, method, pp)].append(filepath)

    for y in ["wA", "w"]:
        order = [pop for pop, sp in sorted(pop2sp.items(), key=lambda kv: sp_sorted(*kv))]
        for (model, sfs, level, method, pp), files in groups.items():
            dico_ada, dico_nn = defaultdict(list), defaultdict(list)
            phylo_ada, phylo_nn = defaultdict(list), defaultdict(list)
            for filepath in files:
                df = pd.read_csv(filepath, sep="\t")
                ddf_ada = df[(df["ADAPTIVE"]) & (df["OMEGA_A"] != "None")]
                ddf_nn = df[(~df["ADAPTIVE"]) & (df["OMEGA_A"] != "None")]

                if len(ddf_nn) == 0 or len(ddf_ada) == 0:
                    continue

                df_phylo = pd.read_csv(filepath.replace(".tsv", ".phylo.tsv"), sep="\t")
                ddf_phylo_ada = df_phylo[df_phylo["ADAPTIVE"]]
                ddf_phylo_nn = df_phylo[~df_phylo["ADAPTIVE"]]

                if y == "wA":
                    dico_ada[y].extend(ddf_ada["OMEGA_A"])
                    dico_nn[y].extend(ddf_nn["OMEGA_A"])
                    phylo_ada[y].extend(ddf_phylo_ada["OMEGA_A"])
                    phylo_nn[y].extend(ddf_phylo_nn["OMEGA_A"])
                else:
                    omega_nn = ddf_nn["OMEGA_A"] + ddf_nn["OMEGA_NA"]
                    ddf_nn = ddf_nn[(0.99 > omega_nn) & (omega_nn > 0.01)]
                    dico_ada[y].extend(ddf_ada["OMEGA_A"] + ddf_ada["OMEGA_NA"])
                    dico_nn[y].extend(ddf_nn["OMEGA_A"] + ddf_nn["OMEGA_NA"])
                    phylo_ada[y].extend(ddf_phylo_ada["OMEGA"])
                    phylo_nn[y].extend(ddf_phylo_nn["OMEGA"])

                _, pop, _, _, _ = filepath.split("/")[-2].replace("_", " ").split("-")
                dico_ada["pop"].extend([format_pop(pop)] * len(ddf_ada))
                dico_nn["pop"].extend([format_pop(pop)] * len(ddf_nn))
                phylo_ada["pop"].extend([format_pop(pop)] * len(ddf_phylo_ada))
                phylo_nn["pop"].extend([format_pop(pop)] * len(ddf_phylo_nn))
            df_ada = pd.DataFrame(dico_ada)
            df_nn = pd.DataFrame(dico_nn)
            df_phylo_ada = pd.DataFrame(phylo_ada)
            df_phylo_nn = pd.DataFrame(phylo_nn)
            sub_pop = set(df_nn["pop"])
            nbr_pop = len(sub_pop)
            sub_order = [pop for pop in order if pop in sub_pop]

            fig, ax = plt.subplots(figsize=(180 * (nbr_pop + 1) / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            graph = sns.violinplot(x="pop", y=y, data=df_nn, ax=ax, color=GREEN, inner="quartile", order=sub_order,
                                   label=f"Nearly-neutral {level}s (subsampled)")
            sns.stripplot(x="pop", y=y, data=df_ada, ax=ax, color=RED, edgecolor="gray", size=10, jitter=0,
                          order=sub_order, label=f"Adaptive {level}s")
            '''
            graph = sns.violinplot(x="pop", y=y, data=df_phylo_nn, ax=ax, color='black', inner="quartile",
                                   order=sub_order)
            sns.stripplot(x="pop", y=y, data=df_phylo_ada, ax=ax, color='black', edgecolor="gray", size=10, jitter=0,
                          order=sub_order)
            '''
            graph.axhline(np.mean(df_phylo_nn[y]), color=GREEN, linestyle="--")
            graph.axhline(np.mean(df_phylo_ada[y]), color=RED, linestyle="--")

            # plt.legend(fontsize=fontsize)

            ax.set_xlabel("", fontsize=fontsize_legend)
            if y == "wA":
                ax.set_ylabel("$\\omega_{\\mathrm{A}}$", fontsize=fontsize)
            else:
                ax.set_ylabel("$d_{\\mathrm{N}} / d_{\\mathrm{S}}$", fontsize=fontsize)
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
            for tick_label in graph.axes.get_xticklabels():
                tick_label.set_fontsize(fontsize_legend)
            for tick_label in graph.axes.get_yticklabels():
                tick_label.set_fontsize(fontsize_legend)
            plt.tight_layout()
            plt.savefig(f"{args.output}/{level}-{method}-{pp}-{sfs}-{model}-{y}.pdf", format="pdf")
            plt.clf()
            plt.close("all")
