import argparse
from glob import glob
import pandas as pd
import numpy as np
from Bio.Phylo.PAML import yn00
from collections import defaultdict
from scipy.optimize import curve_fit


def exp_curve(x, a, b, c):
    return a - b * np.exp(-c * x)


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


omega_dict = defaultdict(list)
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
            yn00_results = read_yn(filepath.replace("_grapes.csv", "_yn00.out"))
            res = list(list(yn00_results.values())[0].values())[0]['YN00']
            dnds = res["omega"]
            omega_dict["OMEGA_A"].append(omega_a)
            omega_dict["ALPHA"].append(omega_a / dnds)
            omega_dict["OMEGA_NA"].append(dnds - omega_a)
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
    elif args.model in ["MK", "aMK"]:
        for filepath in glob(args.folder + "/*.dofe"):
            dofe = open(filepath, 'r')
            dofe.readline()
            dofeline = dofe.readline()
            if "#unfolded" in dofeline:
                dofeline = dofe.readline()
            split_line = dofeline.strip().split("\t")
            Ldn, dn, Lds, ds = split_line[-4:]
            dnds = (int(dn) / float(Ldn)) / (int(ds) / float(Lds))
            nbr_cat = (len(split_line) - 8) // 2
            Lpn, Lps = split_line[2], split_line[nbr_cat + 3]
            assert (Lpn == Ldn)
            assert (Lps == Lds)
            shift = 1
            sfs_n = [int(i) for i in split_line[shift + 3:nbr_cat + 3]]
            sfs_s = [int(i) for i in split_line[nbr_cat + shift + 4:nbr_cat * 2 + 4]]
            if args.model == "MK":
                pnps = (sum(sfs_n) / float(Lpn)) / (sum(sfs_s) / float(Lps))
                omega_dict["P_S"].append(sum(sfs_s) / float(Lps))
                omega_dict["OMEGA_A"].append(dnds - pnps)
                omega_dict["ALPHA"].append((dnds - pnps) / dnds)
                omega_dict["OMEGA_NA"].append(pnps)
                omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
            else:
                assert args.model == "aMK"
                alpha_array = np.array(
                    [((dnds - (pn_x / float(Lpn)) / (ps_x / float(Lps))) / dnds) if ps_x != 0 else np.nan for
                     pn_x, ps_x in zip(sfs_n, sfs_s)])
                mask = np.isfinite(alpha_array)
                delta = 1.0 / len(alpha_array)
                freq_array = np.linspace(delta, 1.0 - delta, len(alpha_array))
                if len(alpha_array[mask]) < 4:
                    continue
                popt, pcov = curve_fit(exp_curve, xdata=freq_array[mask], ydata=alpha_array[mask], check_finite=False,
                                       maxfev=32000, bounds=(0, np.inf))
                aMK = exp_curve(1.0, *popt)
                if aMK < 0.0:
                    continue
                omega_dict["OMEGA_A"].append(aMK * dnds)
                omega_dict["ALPHA"].append(aMK)
                omega_dict["OMEGA_NA"].append((1.0 - aMK) * dnds)
                omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)
    else:
        '''
        from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

        string = ''.join(open("postprocessing.R", "r").readlines())
        postprocessing = SignatureTranslatedAnonymousPackage(string, "postprocessing")
        '''
        for filepath in glob(args.folder + "/*polyDFE.out"):
            polydfe_dico = read_polyDFE(filepath)
            if len(polydfe_dico) == 0: continue
            if 'D_sel' in polydfe_dico and 'D_neut' in polydfe_dico:
                dnds = polydfe_dico["D_sel"] / polydfe_dico["D_neut"]
            else:
                yn00_results = read_yn(filepath.replace("_polyDFE.out", "_yn00.out"))
                res = list(list(yn00_results.values())[0].values())[0]['YN00']
                dnds = res["omega"]
            '''
            estimates = postprocessing.parseOutput(filepath)[0]
            alpha = postprocessing.estimateAlpha(estimates, supLimit=5)[0]
            alpha_dfe = polydfe_dico["alpha_dfe"]
            '''
            alpha = polydfe_dico["alpha_div"]
            omega_dict["OMEGA_NA"].append(dnds * (1 - alpha))
            omega_dict["OMEGA_A"].append(dnds * alpha)
            omega_dict["ALPHA"].append(alpha)
            omega_dict["ADAPTIVE"].append("ADAPTIVE" in filepath)

    pd.DataFrame(omega_dict).to_csv(args.output, sep="\t", index=False)
