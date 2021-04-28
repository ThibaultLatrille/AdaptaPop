import argparse
from libraries import *
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-x', '--xml', required=True, type=str,
                        dest="xml", metavar="<xml>",
                        help="The xml folder")
    parser.add_argument('-s', '--species', required=True, type=str,
                        dest="species", metavar="<species>",
                        help="The species name")

    args = parser.parse_args()

    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, gene_level=True)

    dict_tr_id = build_dict_trID(args.xml, args.species)
    dict_ensg = {v: k for k, v in dict_tr_id.items()}

    ensg_list = sorted(set(dico_omega) & set(dict_ensg))
    df = pd.DataFrame({"ENSG": [k for k in ensg_list], "TR_ID": [dict_ensg[k] for k in ensg_list],
                       "ω": [dico_omega[k][1] for k in ensg_list], "ω0": [dico_omega_0[k][1] for k in ensg_list],
                       "Δω": [dico_omega[k][1] - dico_omega_0[k][1] for k in ensg_list]})

    df.to_csv("TableOmega.tsv", sep="\t", index=False)
