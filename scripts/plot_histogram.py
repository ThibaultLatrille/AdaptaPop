import argparse
import os
import pandas as pd
import statsmodels.api as sm
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

GREEN = "#8FB03E"
RED = "#EB6231"

fontsize = 16
fontsize_legend = 14
my_dpi = 256
fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
bins = 30
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output folder")
    args = parser.parse_args()

    dico_output = {"species": [], "pop": [], "sfs": [], "granularity": [], "model": [],
                   "wA_Selected": [], "wNA_Selected": [], "w_Selected": [], "alpha_Selected": [],
                   "wA_Neutral": [], "wNA_Neutral": [], "w_Neutral": [], "alpha_Neutral": [],
                   "p_val": []}

    for filepath in args.tsv:
        df = pd.read_csv(filepath, sep="\t")
        if len(df["OMEGA_A"]) == 0 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) == 0: continue
        df = df[df["OMEGA_NA"] != "None"]

        sp, pop, granu = filepath.split("/")[-2].replace("_", " ").split("-")
        plot, model, sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")
        name = "{0}/{1}/{2}-{3}-{4}/{5}-{6}".format(args.output, sp.replace(" ", "_"), granu, sfs, model,
                                                    pop.replace(" ", "_"), plot)

        os.makedirs("{0}/{1}/{2}-{3}-{4}".format(args.output, sp.replace(" ", "_"), granu, sfs, model), exist_ok=True)
        assert plot == "histogram"
        # p_val = len([1 for x in nearly_neutral if x > adaptive]) / len(nearly_neutral)
        hist, _, _ = plt.hist(df[~df["ADAPTIVE"]]["OMEGA_A"], bins, density=1, facecolor=GREEN,
                              alpha=0.5, label='Nearly-neutral ({0} subsampling)'.format(sum(~df["ADAPTIVE"])))
        y_max = 1.2 * max(hist)
        hist, _, _ = plt.hist(df[df["ADAPTIVE"]]["OMEGA_A"], bins, density=1, facecolor=RED,
                              alpha=0.5, label='Adaptive ({0} bootstrap)'.format(sum(~df["ADAPTIVE"])))
        y_max = max(y_max, 1.2 * max(hist))
        plt.ylim((0, y_max))

        mean_nearly_neutral = np.mean(df[~df["ADAPTIVE"]]["OMEGA_A"])
        plt.plot((mean_nearly_neutral, mean_nearly_neutral), (0, y_max), linewidth=3, color=GREEN,
                 label=(r'$\omega_{A}$' + '={0:3g}'.format(mean_nearly_neutral)))
        mean_adaptive = np.mean(df[df["ADAPTIVE"]]["OMEGA_A"])
        plt.plot((mean_adaptive, mean_adaptive), (0, y_max), linewidth=3, color=RED,
                 label=(r'$\omega_{A}$' + '={0:3g}'.format(mean_adaptive)))

        smwu, p_val = mannwhitneyu(df[df["ADAPTIVE"]]["OMEGA_A"], df[~df["ADAPTIVE"]]["OMEGA_A"],
                                   alternative="greater")
        plt.title('{0} - {1}\n{2} level, {3} SFS, p-value={4:3g}'.format(sp, pop, granu, sfs, p_val),
                  fontsize=fontsize)
        plt.xlabel(r'$\omega_{A}$', fontsize=fontsize)
        plt.ylabel('Density', fontsize=fontsize)
        plt.legend(fontsize=fontsize_legend, loc='upper left')
        plt.xticks(fontsize=fontsize_legend)
        plt.tight_layout()
        plt.savefig(name + ".pdf", format="pdf")
        plt.savefig(name + ".png", format="png")
        plt.close()

        dico_output["wA_Selected"].append(mean_adaptive)
        dico_output["wNA_Selected"].append(np.mean(df[df["ADAPTIVE"]]["OMEGA_NA"]))
        dico_output["w_Selected"].append(np.mean(df[df["ADAPTIVE"]]["OMEGA_NA"] + df[df["ADAPTIVE"]]["OMEGA_A"]))
        dico_output["alpha_Selected"].append(np.mean(df[df["ADAPTIVE"]]["ALPHA"]))
        dico_output["wA_Neutral"].append(mean_nearly_neutral)
        dico_output["wNA_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["OMEGA_NA"]))
        dico_output["w_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["OMEGA_NA"] + df[~df["ADAPTIVE"]]["OMEGA_A"]))
        dico_output["alpha_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["ALPHA"]))
        dico_output["p_val"].append(p_val)
        dico_output["species"].append(sp)
        dico_output["pop"].append(pop)
        dico_output["granularity"].append(granu)
        dico_output["sfs"].append(sfs)
        dico_output["model"].append(model)
    data = pd.DataFrame(dico_output).to_csv(args.output + "/results.tsv", sep="\t", index=False)
