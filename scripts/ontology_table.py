import scipy.stats as st
from libraries import build_divergence_dico, split_outliers, ontology_table, tex_f
import os
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-x', '--xml', required=True, type=str, dest="xml", metavar="<xml>", help="The xml folder")
    parser.add_argument('-e', '--epistasis', required=False, type=str, default='False', dest="epistasis",
                        metavar="<epistasis>", help="Epistasis (otherwise is adaptive)")
    parser.add_argument('-c', '--ci', required=False, type=str, dest="ci", default="0.025", metavar="<ci>",
                        help="The confidence interval")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-g', '--granularity', required=False, type=str, default=True, dest="granularity",
                        help="Gene or site level")

    args = parser.parse_args()

    args.epistasis = args.epistasis.lower() == 'true'
    gene = args.granularity.lower() == "gene"
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, gene_level=gene, ci=args.ci)
    strg_adap_dico, adap_dico, epi_dico, nn_dico, unclassified_dico = split_outliers(dico_omega_0, dico_omega, gene_level=gene)
    focal_dico = epi_dico if args.epistasis else adap_dico
    go_id2cds_list, go_id2name, set_all_go_cds = ontology_table(args.xml)
    cds_focal_set = set(focal_dico) & set_all_go_cds
    cds_control_set = set(nn_dico) & set_all_go_cds
    if not gene:
        cds_focal_set = set([k for k, v in focal_dico.items() if len(v) / len(dico_omega[k]) >= 0.1]) & set_all_go_cds
        cds_control_set = cds_control_set - cds_focal_set
    assert len(cds_control_set & cds_focal_set) == 0
    cds_set = cds_focal_set.union(cds_control_set)

    print("{0} CDS in focal set.".format(len(cds_focal_set)))
    print("{0} CDS in control set".format(len(cds_control_set)))

    table_tex = open(args.output, 'w')
    table_tex.writelines("\\documentclass[USLetter,5pt]{article}\n\\usepackage{adjustbox}\n")
    table_tex.writelines("\\newcommand{\\specialcell}[2][c]{%\n\\begin{tabular}[#1]{@{}c@{}}#2\\end{tabular}}\n")
    table_tex.writelines("\\begin{document}\n")

    p_value_list = []
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
            oddsratio, p_value = st.fisher_exact(obs, alternative='greater')
            p_value_list.append((go_id, p_value, oddsratio, obs))

    nbr_tests = len(p_value_list)
    e_value_list = [p for p in p_value_list if nbr_tests * p[1] < 1]

    heading = ["Gene Ontology", "$n_{\\mathrm{Observed}}$", "$n_{\\mathrm{Expected}}$", "Odds ratio",
               "$p_{\\mathrm{value}}$", "$e_{\\mathrm{value}}$"]

    table_tex.writelines("\\begin{table}[ht]\n\\centering\n\\begin{adjustbox}{width = 1\\textwidth}\n\\small")
    table_tex.writelines("\\begin{tabular}{" + "|c" * len(heading) + "|}\n")
    table_tex.writelines("\\hline\n")
    table_tex.writelines(" & ".join(heading) + "\\\\\n")
    table_tex.writelines("\\hline\n")

    e_value_list.sort(key=lambda x: x[1])

    for go_id, p_value, oddsratio, obs in e_value_list:
        table_tex.writelines(" & ".join(
            [go_id2name[go_id].replace("_", "-"), str(int(obs[0, 0])), tex_f(obs[0, 0] / oddsratio),
             tex_f(oddsratio), tex_f(p_value), tex_f(p_value * nbr_tests)]) + "\\\\\n")

    caption = "{0} tests performed with {1} CDS detected as {2} and {3} as nearly-neutral." \
              "\n".format(nbr_tests, len(cds_focal_set), "epistasis" if args.epistasis else "adaptive",
                          len(cds_control_set))

    table_tex.writelines("\\hline\n")
    table_tex.writelines("\\end{tabular}\n")
    table_tex.writelines("\\end{adjustbox}\n" + "\\caption{" + caption + "}\n\\end{table}\n")

    table_tex.writelines("\\end{document}\n")
    table_tex.close()
    print('Table generated')
    os.system("pdflatex {0}".format(args.output))
    print('Pdf generated')

    print('Ontology computed')
