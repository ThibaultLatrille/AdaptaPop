import pandas as pd
import scipy.stats as st
from libraries import build_divergence_dico, split_outliers, tex_f, adjusted_holm_pval, format_pval
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
          "$p_{\\mathrm{v}}$", "$p_{\\mathrm{v-adjusted}}$"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-x', '--xml', required=True, type=str, dest="xml", metavar="<xml>", help="The xml folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-g', '--level', required=True, type=str, dest="level", help="Gene or site level")
    parser.add_argument('-m', '--method', required=False, type=str, default='MutSel', dest="method",
                        help="Method to detect adaptation (MutSel, Classical, MutSelExclu)")
    parser.add_argument('-p', '--pp', required=True, type=str, dest="pp", help="Posterior probability")
    args = parser.parse_args()

    gene = args.level.lower() == "gene"
    list_ensg = [i[:-3] for i in os.listdir(args.folder)]
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, list_ensg, gene_level=gene, pp=args.pp)
    adap_dico, nn_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=gene, method=args.method)
    go_id2cds_list, go_id2name, set_all_go_cds = ontology_table(args.xml)

    if gene:
        cds_focal_set = set(adap_dico) & set_all_go_cds
        # cds_control_set = set(nn_dico) & set_all_go_cds
    else:
        cds_focal_set = set()
        for ensg, site_list in adap_dico.items():
            if ensg in set_all_go_cds and len(site_list) > 1:
                cds_focal_set.add(ensg)
        # dico_omega_0, dico_omega = build_divergence_dico(args.folder, list_ensg, gene_level=True)
        # _, nn_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=True, method=args.method)
        # cds_control_set = (set(nn_dico) & set_all_go_cds) - cds_focal_set
    cds_control_set = set_all_go_cds - cds_focal_set
    assert len(cds_control_set & cds_focal_set) == 0
    cds_set = cds_focal_set.union(cds_control_set)

    print("{0} {1}s in focal set.".format(len(cds_focal_set), args.level))
    print("{0} {1}s in control set".format(len(cds_control_set), args.level))

    dico_ouput = {"GO": [], "obs": [], "exp": [], "oddsratio": [], "pval": []}
    for go_id, cds_all_go_set in go_id2cds_list.items():
        cds_go_set = cds_all_go_set & cds_set
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
            oddsratio, pval = st.fisher_exact(obs, alternative='greater')
            dico_ouput["GO"].append(go_id2name[go_id].replace("_", "-").replace("[", "").replace("]", ""))
            dico_ouput["obs"].append(int(obs[0, 0]))
            dico_ouput["exp"].append(float(obs[0, 0]) / oddsratio)
            dico_ouput["oddsratio"].append(oddsratio)
            dico_ouput["pval"].append(pval)

    df = pd.DataFrame(dico_ouput)
    df = adjusted_holm_pval(df, alpha=0.05, format_p=False)
    df.sort_values(by="pval_adj", inplace=True)
    df.to_csv(args.output.replace(".tex", ".tsv"), index=False, sep="\t")

    text_core = f"{len(df['pval'])} tests performed with {len(cds_focal_set)} genes detected with {args.method.lower()}, and {len(cds_control_set)} as control."

    text_core += "\\scriptsize\n"
    df_head = format_pval(df[df["pval_adj"] < 1.0]).head(500)
    text_core += df_head.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=header,
                                  column_format="|l|r|r|r|r|r|")

    with open(args.output.replace(".tex", "") + ".core.tex", 'w') as core_f:
        core_f.write(text_core)

    with open(args.output, "w") as table_tex:
        table_tex.write("\n".join(open(os.path.dirname(args.output) + "/table-import.tex", 'r').readlines()))
        table_tex.write("\\begin{document}\n")
        table_tex.write(text_core)
        table_tex.write("\\end{document}\n")

    print('Table generated')
    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {1}".format(
        os.path.dirname(args.output),
        args.output)
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
    print('Pdf generated')

    print('Ontology computed')
