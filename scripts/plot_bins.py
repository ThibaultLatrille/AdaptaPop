import argparse
from glob import glob
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

GREEN = "#8FB03E"

my_dpi = 256
fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
bins = 30
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    omega_dict = {"OMEGA_NA": [], "OMEGA_A": [], "OMEGA_0": [], "ALPHA": []}
    for filepath in sorted(glob(args.folder + "/*.csv")):
        if args.model in ["dfem", "grapes"]:
            dfem_df = pd.read_csv(filepath)
            omega_a = float(dfem_df[dfem_df["model"] == "GammaExpo"]["omegaA"])
            alpha = float(dfem_df[dfem_df["model"] == "GammaExpo"]["alpha"])
            omega_dict["OMEGA_A"].append(omega_a)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["OMEGA_NA"].append(omega_a * (1 - alpha) / alpha)
            omega_dict["OMEGA_0"].append(float(filepath.split("/")[-1].split("_")[0]))

    pd.DataFrame(omega_dict).to_csv(args.output.replace(".pdf", ".tsv"), sep="\t", index=False, float_format="%.3g")
    plt.scatter(omega_dict["OMEGA_0"], omega_dict["OMEGA_NA"])

    idf = np.linspace(min(omega_dict["OMEGA_0"]) - 0.05, max(omega_dict["OMEGA_0"]) + 0.05, 30)
    model = sm.OLS(list(omega_dict["OMEGA_NA"]), sm.add_constant(list(omega_dict["OMEGA_0"])))
    results = model.fit()
    b, a = results.params[0:2]
    plt.plot(idf, a * idf + b, 'r-',
             label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(float(a), abs(float(b)), results.rsquared,
                                                                      "+" if float(b) > 0 else "-"), color=GREEN)

    plt.xlabel(r'$\omega_{NA}$', fontsize=16)
    plt.ylabel(r'$\omega_{0}$', fontsize=16)
    plt.tight_layout()
    plt.legend(fontsize=14, loc="upper left")
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.close()
