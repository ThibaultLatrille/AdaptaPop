import argparse
from glob import glob
import pandas as pd

omega_dict = {"OMEGA_NA": [], "OMEGA_A": [], "OMEGA_0": [], "ALPHA": [], "ADAPTIVE": []}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    for filepath in glob(args.folder + "/*.csv"):
        if args.model in ["dfem", "grapes"]:
            dfem_df = pd.read_csv(filepath)
            omega_a = float(dfem_df[dfem_df["model"] == "GammaExpo"]["omegaA"])
            alpha = float(dfem_df[dfem_df["model"] == "GammaExpo"]["alpha"])
            omega_dict["OMEGA_A"].append(omega_a)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["OMEGA_NA"].append(omega_a * (1 - alpha) / alpha if alpha != 0 else "NaN")
            omega_dict["OMEGA_0"].append(float(filepath.split("/")[-1].split("_")[0]) if "Bins" in filepath else "NaN")
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)

    pd.DataFrame(omega_dict).to_csv(args.output, sep="\t", index=False)
