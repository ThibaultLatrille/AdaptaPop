import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output folder")
    args = parser.parse_args()

    dico_output = defaultdict(list)
    dico_ps = defaultdict(list)

    for filepath in args.tsv:
        df = pd.read_csv(filepath, sep="\t")
        if len(df["OMEGA_A"]) == 0 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) == 0: continue
        df = df[df["OMEGA_NA"] != "None"]

        sp, pop, granu = filepath.split("/")[-2].replace("_", " ").split("-")
        plot, model, sfs = filepath.split("/")[-1].replace(".tsv", "").split("-")
        wa_ada = np.mean(df[df["ADAPTIVE"]]["OMEGA_A"])
        p_val = len([1 for x in df[~df["ADAPTIVE"]]["OMEGA_A"] if x > wa_ada]) / len(df[~df["ADAPTIVE"]]["OMEGA_A"])
        wa_nn_mean = np.mean(df[~df["ADAPTIVE"]]["OMEGA_A"])
        dico_output["wA_Selected"].append(wa_ada)
        dico_output["wA_Neutral"].append(wa_nn_mean)
        dico_output["DeltawA_Neutral"].append(wa_ada - wa_nn_mean)
        # dico_output["wNA_Selected"].append(np.mean(df[df["ADAPTIVE"]]["OMEGA_NA"]))
        # dico_output["w_Selected"].append(np.mean(df[df["ADAPTIVE"]]["OMEGA_NA"] + df[df["ADAPTIVE"]]["OMEGA_A"]))
        # dico_output["alpha_Selected"].append(np.mean(df[df["ADAPTIVE"]]["ALPHA"]))
        # dico_output["wNA_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["OMEGA_NA"]))
        # dico_output["w_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["OMEGA_NA"] + df[~df["ADAPTIVE"]]["OMEGA_A"]))
        # dico_output["alpha_Neutral"].append(np.mean(df[~df["ADAPTIVE"]]["ALPHA"]))
        dico_output["p_val"].append(p_val)
        dico_output["species"].append(sp)
        dico_output["pop"].append(pop)
        dico_output["granularity"].append(granu)
        dico_output["sfs"].append(sfs)
        dico_output["model"].append(model)
        if "P_S" in df:
            dico_ps["P_S"].append(np.mean(df["P_S"]))
            dico_ps["pop"].append(pop)
    data = pd.DataFrame(dico_output)
    dfps = pd.DataFrame(dico_ps).groupby("pop").mean()

    merge = pd.merge(data, dfps, how="inner", on=["pop"])
    merge.to_csv(args.output + "/results.tsv", sep="\t", index=False)
