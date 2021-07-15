import pandas as pd
import scipy.stats as st
from libraries import build_divergence_dico, split_outliers, tex_f
import os
from lxml import etree
import numpy as np
import argparse


def ontology_table(xml_folder):
    print('Finding CDS ontologies.')
    go_id2name, go_id2cds_list = {}, {}
    all_go_set = set()
    for file in os.listdir(xml_folder):
        root = etree.parse(xml_folder + "/" + file).getroot()
        for annot in root.find('goAnnot').findall("annot"):
            go_id = annot.find('goId').text
            go_name = annot.find('goName').text
            if go_id not in go_id2name: go_id2name[go_id] = go_name.replace('"', '')
            if go_id not in go_id2cds_list: go_id2cds_list[go_id] = set()
            ensg = file.replace(".xml", "")
            go_id2cds_list[go_id].add(ensg)
            all_go_set.add(ensg)
    print('CDS ontologies found.')
    return go_id2cds_list, go_id2name, all_go_set


header = ["Gene Ontology", "$n_{\\mathrm{Observed}}$", "$n_{\\mathrm{Expected}}$", "Odds ratio",
          "$p_{\\mathrm{value}}$", "$e_{\\mathrm{value}}$"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-x', '--xml', required=True, type=str, dest="xml", metavar="<xml>", help="The xml folder")
    parser.add_argument('-c', '--category', required=False, type=str, default='adaptive', dest="category",
                        help="Category of genes or sites (either epistasis, adaptive or strongly_adaptive)")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-g', '--granularity', required=True, type=str, default="gene", dest="granularity",
                        help="Gene or site level")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    list_ensg = [i[:-3] for i in os.listdir(args.folder)]
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, list_ensg, gene_level=gene)
    strg_adap_dico, adap_dico, epi_dico, nn_dico, unclassified_dico = split_outliers(dico_omega_0, dico_omega,
                                                                                     gene_level=gene)
    focal_dico = epi_dico if args.category.lower() == "epistasis" else (
        adap_dico if args.category.lower() == "adaptive" else strg_adap_dico)
    go_id2cds_list, go_id2name, set_all_go_cds = ontology_table(args.xml)
    if gene:
        cds_focal_set = set(focal_dico) & set_all_go_cds
        cds_control_set = set(nn_dico) & set_all_go_cds
    else:
        cds_focal_set, cds_control_set = set(), set()
        for k, v in focal_dico.items():
            if k not in set_all_go_cds: continue
            cds_focal_set.update(set([k + str(i) for i in v]))
        for k, v in nn_dico.items():
            if k not in set_all_go_cds: continue
            cds_control_set.update(set([k + str(i) for i in v]))
    assert len(cds_control_set & cds_focal_set) == 0
    cds_set = cds_focal_set.union(cds_control_set)

    print("{0} {1}s in focal set.".format(len(cds_focal_set), args.granularity))
    print("{0} {1}s in control set".format(len(cds_control_set), args.granularity))

    dico_ouput = {"GO": [], "obs": [], "exp": [], "oddsratio": [], "p_value": []}
    for go_id, cds_all_go_set in go_id2cds_list.items():
        if gene:
            cds_go_set = cds_all_go_set & cds_set
        else:
            cds_go_set = set()
            for cds_go in cds_all_go_set:
                if cds_go in focal_dico:
                    cds_go_set.update(set([cds_go + str(i) for i in focal_dico[cds_go]]))
                if cds_go in nn_dico:
                    cds_go_set.update(set([cds_go + str(i) for i in nn_dico[cds_go]]))
            print(go_id2name[go_id].replace("_", "-"))
        cds_no_go_set = cds_set - cds_go_set
        obs = np.array([[len(cds_go_set & cds_focal_set),
                         len(cds_no_go_set & cds_focal_set)],
                        [len(cds_go_set & cds_control_set),
                         len(cds_no_go_set & cds_control_set)]])
        assert sum(obs[0, :]) == len(cds_focal_set)
        assert sum(obs[1, :]) == len(cds_control_set)
        assert sum(obs[:, 0]) == len(cds_go_set)
        assert sum(obs[:, 1]) == len(cds_no_go_set)
        assert np.sum(obs) == len(cds_set)

        if np.min(obs) > 1:
            oddsratio, p_value = st.fisher_exact(obs, alternative='greater')
            dico_ouput["GO"].append(go_id2name[go_id].replace("_", "-").replace("[acyl-carrier-protein]", ""))
            dico_ouput["obs"].append(int(obs[0, 0]))
            dico_ouput["exp"].append(float(obs[0, 0]) / oddsratio)
            dico_ouput["oddsratio"].append(oddsratio)
            dico_ouput["p_value"].append(p_value)

    df = pd.DataFrame(dico_ouput)
    df["e_value"] = df["p_value"] * len(df["p_value"])
    df.sort_values(by="e_value", inplace=True)
    df.to_csv(args.output.replace(".tex", ".tsv"), index=False, sep="\t")

    table_tex = open(args.output, 'w')
    latex_in = open(os.path.dirname(args.output) + "/table-import.tex", 'r')
    table_tex.write("\n".join(latex_in.readlines()))
    latex_in.close()

    table_tex.write("\\begin{document}\n")
    table_tex.write("{0} tests performed with {1} {4}s detected as {2} and {3} as nearly-neutral." \
                    "\n".format(len(df["p_value"]), len(cds_focal_set), args.category.lower().replace("-", " "),
                                len(cds_control_set), args.granularity))
    table_tex.write("\\scriptsize\n")
    table_tex.write(df.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=header))
    table_tex.write("\\end{document}\n")
    table_tex.close()
    print('Table generated')
    os.system("pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {1}".format(os.path.dirname(args.output), args.output))
    print('Pdf generated')

    print('Ontology computed')
