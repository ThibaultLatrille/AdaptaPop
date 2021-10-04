import argparse
from glob import glob
import pandas as pd
from Bio.Phylo.PAML import yn00


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


def read_yn(path):
    return yn00.read(path)


omega_dict = {"OMEGA_NA": [], "OMEGA_A": [], "ALPHA": [], "ADAPTIVE": []}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    if args.model in ["dfem", "grapes"]:
        for filepath in glob(args.folder + "/*{0}.csv".format(args.model)):
            dfem_df = pd.read_csv(filepath)
            ge_df = dfem_df[dfem_df["model"] == "GammaExpo"]
            omega_a = float(ge_df["omegaA"])
            alpha = float(ge_df["alpha"])
            omega_dict["OMEGA_A"].append(omega_a)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["OMEGA_NA"].append(omega_a * (1 - alpha) / alpha if alpha != 0 else "NaN")
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
    elif args.model == "MK":
        # for filepath in glob(args.folder + "/*grapes.csv"):
        #     dfem_df = pd.read_csv(filepath)
        #     df = dfem_df[dfem_df["model"] == "GammaExpo"]
        for filepath in glob(args.folder + "/*{0}.tsv".format(args.model)):
            df = pd.read_csv(filepath, sep='\t')
            dnds = (float(df["dn"]) / float(df["Ldn"])) / (float(df["ds"]) / float(df["Lds"]))
            pnps = (float(df["pn"]) / float(df["Lpn"])) / (float(df["ps"]) / float(df["Lps"]))
            omega_dict["OMEGA_A"].append(dnds - pnps)
            omega_dict["ALPHA"].append((dnds - pnps) / dnds)
            omega_dict["OMEGA_NA"].append(pnps)
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
    else:
        from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

        string = ''.join(open("postprocessing.R", "r").readlines())
        postprocessing = SignatureTranslatedAnonymousPackage(string, "postprocessing")

        for filepath in glob(args.folder + "/*polyDFE.out"):
            polydfe_dico = read_polyDFE(filepath)
            if len(polydfe_dico) == 0: continue
            if 'D_sel' in polydfe_dico and 'D_neut' in polydfe_dico:
                omega = polydfe_dico["D_sel"] / polydfe_dico["D_neut"]
            else:
                yn00_results = read_yn(filepath.replace("_polyDFE.out", "_yn00.out"))
                res = list(list(yn00_results.values())[0].values())[0]['YN00']
                omega = res["omega"]

            # estimates = postprocessing.parseOutput(filepath)[0]
            # alpha = postprocessing.estimateAlpha(estimates, supLimit=5)[0]
            # alpha_dfe = polydfe_dico["alpha_dfe"]
            alpha = polydfe_dico["alpha_div"]
            omega_dict["OMEGA_NA"].append(omega * (1 - alpha))
            omega_dict["OMEGA_A"].append(omega * alpha)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)

    pd.DataFrame(omega_dict).to_csv(args.output, sep="\t", index=False)
