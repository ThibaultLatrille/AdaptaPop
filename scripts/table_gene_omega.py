import argparse
import os
from libraries import build_divergence_dico, split_outliers, build_dict_trID
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-x', '--xml', required=True, type=str, dest="xml", help="The xml folder")
    parser.add_argument('-s', '--species', required=True, type=str, nargs="+", dest="species", help="The species name")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()

    ensg_list = sorted([i[:-3] for i in os.listdir(args.div_folder) if i.startswith("ENSG")])
    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, ensg_list, gene_level=True)
    adap, nn, _ = split_outliers(dico_omega_0, dico_omega, gene_level=True)
    dico_cat = dict()
    for ensg_cat_list, cat in [(adap, "adaptive"), (nn, "nearly-neutral")]:
        for ensg in ensg_cat_list:
            dico_cat[ensg] = cat

    dico_out = {"ENSG": [k for k in ensg_list],
                "ω_lower": [dico_omega[k][0] for k in ensg_list],
                "ω": [dico_omega[k][1] for k in ensg_list],
                "ω_upper": [dico_omega[k][2] for k in ensg_list],
                "ω0_lower": [dico_omega_0[k][0] for k in ensg_list],
                "ω0": [dico_omega_0[k][1] for k in ensg_list],
                "ω0_upper": [dico_omega_0[k][2] for k in ensg_list],
                "ωA_phy": [dico_omega[k][1] - dico_omega_0[k][1] for k in ensg_list],
                "category": [dico_cat[k] if k in dico_cat else "unclassified" for k in ensg_list]}

    for sp in args.species:
        dict_tr_id = build_dict_trID(args.xml, sp)
        dict_ensg = {v: k for k, v in dict_tr_id.items()}
        dico_out["TR_ID_" + sp] = [dict_ensg[k] if k in dict_ensg else "Null" for k in ensg_list]

    pd.DataFrame(dico_out).to_csv(args.output, index=False, sep=",")
