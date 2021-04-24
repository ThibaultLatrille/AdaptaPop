import argparse
from glob import glob
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

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

    omega_dict = {"OMEGA_A": [], "BIN": []}
    for filepath in glob(args.folder + "/*.csv"):

        dfem_df = pd.read_csv(filepath)
        omega_a = float(dfem_df[dfem_df["model"] == args.model]["omegaA"])

        omega_dict["OMEGA_A"].append(omega_a)
        omega_dict["BIN"].append(filepath.split("/")[-1].split("_")[0])

    pd.DataFrame(omega_dict).to_csv(args.output.replace(".pdf", ".tsv"), sep="\t", index=False)
    # p_val = len([1 for x in nearly_neutral if x > adaptive]) / len(nearly_neutral)
    plt.scatter(omega_dict["BIN"], omega_dict["OMEGA_A"])
    # plt.title('p-value={0:3g}'.format(p_val), fontsize=8)
    plt.xlabel(r'$\omega_{A}$', fontsize=8)
    plt.ylabel(r'$\omega_{A}$', fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.close()
