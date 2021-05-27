import argparse
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

    dico_output = {"species": [], "pop": [], "sfs": [], "granularity": [],
                   "wA_Adaptive": [], "wA_Nearly_Neutral": [], "p_val": [], "a": [], "r2": []}

    for filepath in args.tsv:
        df = pd.read_csv(filepath, sep="\t")
        if len(df["OMEGA_A"]) < 5 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) < 5: continue
        df = df[df["OMEGA_NA"] != "None"]

        sp, pop, granu = filepath.split("/")[-2].replace("_", " ").split("-")
        sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")[-1]
        name = "{0}/{1}-{2}".format(args.output, filepath.split("/")[-2], filepath.split("/")[-1].replace(".tsv", ""))

        dico_output["species"].append(sp)
        dico_output["pop"].append(pop)
        dico_output["granularity"].append(granu)
        dico_output["sfs"].append(sfs)

        if "histogram" in filepath:
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
            plt.savefig(name + "-hist.pdf", format="pdf")
            plt.savefig(name + "-hist.pdf", format="png")
            plt.close()

            dico_output["wA_Adaptive"].append(mean_adaptive)
            dico_output["wA_Nearly_Neutral"].append(mean_nearly_neutral)
            dico_output["p_val"].append(p_val)
            dico_output["a"].append("NaN")
            dico_output["r2"].append("NaN")

        elif "bins" in filepath:
            plt.scatter(df["OMEGA_0"], df["OMEGA_NA"])

            idf = np.linspace(min(df["OMEGA_0"]) - 0.05, max(df["OMEGA_0"]) + 0.05, 30)
            model = sm.OLS(list(df["OMEGA_NA"]), sm.add_constant(list(df["OMEGA_0"])))
            results = model.fit()
            b, a = results.params[0:2]
            plt.plot(idf, a * idf + b, 'r-',
                     label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(float(a), abs(float(b)), results.rsquared,
                                                                              "+" if float(b) > 0 else "-"),
                     color=GREEN)

            plt.xlabel(r'$\omega_{NA}$', fontsize=16)
            plt.ylabel(r'$\omega_{0}$', fontsize=16)
            plt.title('{0} - {1}\n{2} level, {3} SFS'.format(sp, pop, granu, sfs), fontsize=fontsize)
            plt.tight_layout()
            plt.legend(fontsize=14, loc="upper left")
            plt.savefig(name + "-scatter.pdf", format="pdf")
            plt.savefig(name + "-scatter.png", format="png")
            plt.close()

            dico_output["wA_Adaptive"].append("NaN")
            dico_output["wA_Nearly_Neutral"].append("NaN")
            dico_output["p_val"].append("NaN")
            dico_output["a"].append(a)
            dico_output["r2"].append(results.rsquared)
    data = pd.DataFrame(dico_output).to_csv(args.output + "/results.tsv", sep="\t", index=False)
