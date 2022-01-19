import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from libraries import sort_df, tex_f

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, nargs="+", dest="tsv", help="Input tsv files")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tsv file")
    parser.add_argument('-s', '--sample_list', required=False, type=str, dest="sample_list", help="Sample list file")
    args = parser.parse_args()

    dico_output = defaultdict(list)
    dico_ps = defaultdict(list)

    for filepath in args.tsv:
        df = pd.read_csv(filepath, sep="\t")
        if len(df["OMEGA_A"]) == 0 or len(df[~df["ADAPTIVE"]]["OMEGA_A"]) == 0:
            continue
        df_phylo = pd.read_csv(filepath.replace(".tsv", ".phylo.tsv"), sep="\t")
        df = df[df["OMEGA_NA"] != "None"]

        sp, pop, model, sfs = filepath.split("/")[-1].replace("_", " ").replace(".tsv", "").split("-")

        dico_output["species"].append(sp)
        dico_output["pop"].append(pop)
        dico_output["sfs"].append(sfs)
        dico_output["model"].append(model)
        dico_output["ωA"].append(np.mean(df["OMEGA_A"]))
        dico_output["Δω"].append(np.mean(df_phylo["OMEGA_A"]))

    data = pd.DataFrame(dico_output)
    data = sort_df(data, args.sample_list)
    data.to_csv(args.output, sep="\t", index=False)
    data.to_latex(args.output.replace(".tsv", ".tex"), index=False, escape=False, float_format=tex_f,
                  column_format="|l|l|r|r|r|r|")
