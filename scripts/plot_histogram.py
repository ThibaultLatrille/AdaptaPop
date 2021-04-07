import argparse
from glob import glob
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

GREEN = "#8FB03E"
RED = "#EB6231"

nearly_neutral = []
adaptive = []
omega_dict = {"OMEGA_A": [], "ADAPTIVE": []}

my_dpi = 256
fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
bins = 30
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    for filepath in glob(args.folder + "/*.csv"):

        dfem_df = pd.read_csv(filepath)
        omega_a = float(dfem_df[dfem_df["model"] == args.model]["omegaA"])

        omega_dict["OMEGA_A"].append(omega_a)

        if "ADAPTIVE" in filepath:
            omega_dict["ADAPTIVE"].append(True)
            adaptive.append(omega_a)
        else:
            omega_dict["ADAPTIVE"].append(False)
            nearly_neutral.append(omega_a)

    pd.DataFrame(omega_dict).to_csv(args.output.replace(".pdf", ".tsv"), sep="\t", index=False)
    # p_val = len([1 for x in nearly_neutral if x > adaptive]) / len(nearly_neutral)
    hist, _, _ = plt.hist(nearly_neutral, bins, density=1, facecolor=GREEN,
                          alpha=0.4, label='Nearly-neutral ({0} resampling)'.format(len(nearly_neutral)))
    y_max = 1.2 * max(hist)
    hist, _, _ = plt.hist(adaptive, bins, density=1, facecolor=RED,
                          alpha=0.4, label='Adaptive ({0} resampling)'.format(len(adaptive)))
    y_max = max(y_max, 1.2 * max(hist))
    plt.ylim((0, y_max))
    # plt.title('p-value={0:3g}'.format(p_val), fontsize=8)
    mean_nearly_neutral = np.mean(nearly_neutral)
    plt.plot((mean_nearly_neutral, mean_nearly_neutral), (0, y_max), linewidth=3, color=GREEN,
             label=(r'$\omega_{A}$' + '={0:3g}'.format(mean_nearly_neutral)))
    mean_adaptive = np.mean(adaptive)
    plt.plot((mean_adaptive, mean_adaptive), (0, y_max), linewidth=3, color=RED,
             label=(r'$\omega_{A}$' + '={0:3g}'.format(mean_adaptive)))
    plt.xlabel(r'$\omega_{A}$', fontsize=8)
    plt.ylabel('Density', fontsize=8)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.close()
