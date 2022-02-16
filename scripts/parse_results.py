import argparse
from glob import glob
import pandas as pd
import numpy as np
from Bio.Phylo.PAML import yn00
from collections import defaultdict
from libraries import theta


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


def parse_line(line):
    return float(line.split(",")[0].split("=")[-1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=False, type=str, dest="folder", help="Folder path")
    parser.add_argument('-m', '--model', required=False, type=str, dest="model", help="Model name")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    omega_dict = defaultdict(list)
    phylo_dict = defaultdict(list)
    for filepath in glob(args.folder + "/*.txt"):
        f = open(filepath, 'r')
        while True:
            omega_line = f.readline()
            if omega_line.strip() == "":
                break
            omega_0_line = f.readline()
            omega_A_line = f.readline()
            if "all adaptive" in omega_line:
                adaptive = True
            elif ("resampled nearly" in omega_line) or ('all genome' in omega_line):
                adaptive = False
            else:
                continue
            phylo_dict["OMEGA"].append(parse_line(omega_line))
            phylo_dict["OMEGA_0"].append(parse_line(omega_0_line))
            phylo_dict["OMEGA_A"].append(parse_line(omega_A_line))
            phylo_dict["ADAPTIVE"].append(adaptive)

    pd.DataFrame(phylo_dict).to_csv(args.output.replace(".tsv", ".phylo.tsv"), sep="\t", index=False)

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
    elif args.model in ["MK", "MK_tajima", "MK_fay_wu"]:
        for filepath in glob(args.folder + "/*.dofe"):
            dofe = open(filepath, 'r')
            dofe.readline()
            dofeline = dofe.readline()
            if "#unfolded" in dofeline:
                dofeline = dofe.readline()
            split_line = dofeline.strip().split("\t")
            Ldn, dn, Lds, ds = [float(i) for i in split_line[-4:]]
            dnds = (dn / Ldn) / (ds / Lds)
            nbr_cat = (len(split_line) - 8) // 2
            Lpn, Lps = float(split_line[2]), float(split_line[nbr_cat + 3])
            assert (Lpn == Ldn)
            assert (Lps == Lds)
            shift = 1
            sfs_n = [int(i) for i in split_line[shift + 3:nbr_cat + 3]]
            sfs_s = [int(i) for i in split_line[nbr_cat + shift + 4:nbr_cat * 2 + 4]]
            if args.model == "MK":
                ps = sum(sfs_s) / Lps
                pn = sum(sfs_n) / Lpn
                pnps = pn / ps
                ps = theta(sfs_s, nbr_cat, "watterson") / Lps
                pn = theta(sfs_n, nbr_cat, "watterson") / Lpn
                assert abs(pn / ps - pnps) < 1e-6
            elif args.model == "MK_tajima":
                ps = theta(sfs_s, nbr_cat, "tajima") / Lps
                pn = theta(sfs_n, nbr_cat, "tajima") / Lpn
            else:
                assert args.model == "MK_fay_wu"
                ps = theta(sfs_s, nbr_cat, "fay_wu") / Lps
                pn = theta(sfs_n, nbr_cat, "fay_wu") / Lpn
            pnps = pn / ps
            omega_dict["Ldn"].append(Ldn)
            omega_dict["dn"].append(dn)
            omega_dict["Lds"].append(Lds)
            omega_dict["ds"].append(ds)
            omega_dict["P_S"].append(ps)
            omega_dict["OMEGA_A"].append(dnds - pnps)
            omega_dict["ALPHA"].append((dnds - pnps) / dnds)
            omega_dict["OMEGA_NA"].append(pnps)
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
