import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tsv file")
    args = parser.parse_args()

    dico_output = defaultdict(list)
    dico_ps = defaultdict(list)

    for filepath in args.tsv:
        df = pd.read_csv(filepath, sep="\t")
        if len(df["OMEGA_A"]) == 0 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) == 0:
            continue
        df_phylo = pd.read_csv(filepath.replace(".tsv", ".phylo.tsv"), sep="\t")
        df = df[df["OMEGA_NA"] != "None"]

        sp, pop, level, method, pp = filepath.split("/")[-2].replace("_", " ").split("-")[:5]
        plot, model, sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")[:3]

        wa_test = np.mean(df[df["ADAPTIVE"]]["OMEGA_A"])
        wa_pval = len([1 for x in df[~df["ADAPTIVE"]]["OMEGA_A"] if x > wa_test]) / len(df[~df["ADAPTIVE"]]["OMEGA_A"])

        w_test = np.mean(df[df["ADAPTIVE"]]["OMEGA_NA"] + df[df["ADAPTIVE"]]["OMEGA_A"])
        w_control = df[~df["ADAPTIVE"]]["OMEGA_A"] + df[~df["ADAPTIVE"]]["OMEGA_NA"]
        w_pval = len([1 for x in w_control if x > w_test]) / len(w_control)

        wa_control_mean = np.mean(df[~df["ADAPTIVE"]]["OMEGA_A"])

        dico_output["wA_Test"].append(wa_test)
        dico_output["wA_Control"].append(wa_control_mean)
        dico_output["wA_pval"].append(wa_pval)
        dico_output["w_Test"].append(w_test)
        dico_output["w_Control"].append(np.mean(w_control))
        dico_output["w_pval"].append(w_pval)
        dico_output["Δw_Test"].append(np.mean(df_phylo[df_phylo["ADAPTIVE"]]["OMEGA_A"]))
        dico_output["Δw_Control"].append(np.mean(df_phylo[~df_phylo["ADAPTIVE"]]["OMEGA_A"]))
        dico_output["species"].append(sp)
        dico_output["pop"].append(pop)
        dico_output["level"].append(level)
        dico_output["method"].append(method)
        dico_output["pp"].append(pp)
        dico_output["sfs"].append(sfs)
        dico_output["model"].append(model)
        if "P_S" in df:
            dico_ps["P_S"].append(np.mean(df["P_S"]))
            dico_ps["pop"].append(pop)

    data = pd.DataFrame(dico_output)
    dfps = pd.DataFrame(dico_ps).groupby("pop").mean()

    merge = pd.merge(data, dfps, how="inner", on=["pop"])
    merge.to_csv(args.output, sep="\t", index=False)
