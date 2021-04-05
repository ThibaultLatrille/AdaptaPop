import argparse
from glob import glob
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

GREEN = "#8FB03E"
RED = "#EB6231"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    non_outliers = []
    outliers = 0
    omega_dict = {"OMEGA_A": [], "OUTLIERS": [], "REP": []}

    for filepath in glob(args.folder + "/*.csv"):

        dfem_df = pd.read_csv(filepath)
        omega_a = float(dfem_df[dfem_df["model"] == args.model]["omegaA"])

        omega_dict["OMEGA_A"].append(omega_a)

        if "OUTLIERS" in filepath:
            omega_dict["OUTLIERS"].append(True)
            omega_dict["REP"].append(None)
            outliers = omega_a
        else:
            omega_dict["OUTLIERS"].append(False)
            omega_dict["REP"].append(filepath.split("/")[-1].split(".")[0])
            non_outliers.append(omega_a)

    pd.DataFrame(omega_dict).to_csv(args.output.replace(".pdf", ".tsv"), sep="\t", index=False)
    my_dpi = 256
    fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    p_val = len([1 for x in non_outliers if x > outliers]) / len(non_outliers)

    bins = 50
    hist, _, _ = plt.hist(non_outliers, bins, density=1, facecolor=GREEN,
                          alpha=0.4, label='Non-outliers CDS ({0} resampling)'.format(len(non_outliers)))
    y_max = 1.2 * max(hist)
    plt.ylim((0, y_max))
    plt.title('p-value={0:3g}'.format(p_val), fontsize=8)
    plt.plot((outliers, outliers), (0, y_max), linewidth=3, color=RED,
             label=(r'$\omega_{A}$' + '={0:3g})'.format(outliers)))
    plt.xlabel(r'$\omega_{A}$', fontsize=8)
    plt.ylabel('Density', fontsize=8)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.close()
