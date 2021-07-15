import argparse
from glob import glob
import pandas as pd


def read_polyDFE(path):
    dico_out = {}
    with open(path, 'r') as file:
        for line in file.readlines():
            if line.count("=") != 1: continue
            k, v = line.strip().split("=")
            for ch in ["E[", "]", "-", " "]:
                k = k.replace(ch, "")
            dico_out[k] = float(v)
    return dico_out


omega_dict = {"OMEGA_NA": [], "OMEGA_A": [], "OMEGA_0": [], "ALPHA": [], "ADAPTIVE": []}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    if args.model in ["dfem", "grapes"]:
        for filepath in glob(args.folder + "/*.csv"):
            dfem_df = pd.read_csv(filepath)
            omega_a = float(dfem_df[dfem_df["model"] == "GammaExpo"]["omegaA"])
            alpha = float(dfem_df[dfem_df["model"] == "GammaExpo"]["alpha"])
            omega_dict["OMEGA_A"].append(omega_a)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["OMEGA_NA"].append(omega_a * (1 - alpha) / alpha if alpha != 0 else "NaN")
            omega_dict["OMEGA_0"].append(float(filepath.split("/")[-1].split("_")[0]) if "Bins" in filepath else "NaN")
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
    else:
        for filepath in glob(args.folder + "/*polyDFE.out"):
            polydfe_dico = read_polyDFE(filepath)
            omega = polydfe_dico["D_sel"] / polydfe_dico["D_neut"]
            alpha = polydfe_dico["alpha_div"]
            omega_dict["OMEGA_NA"].append(omega * (1 - alpha))
            omega_dict["OMEGA_A"].append(omega * alpha)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["OMEGA_0"].append(float(filepath.split("/")[-1].split("_")[0]) if "Bins" in filepath else "NaN")
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)

    pd.DataFrame(omega_dict).to_csv(args.output, sep="\t", index=False)
