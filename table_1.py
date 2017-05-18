import numpy as np
import os
import scipy.stats as st
from cds_libraries import load_pb_table, str_to_table
import pickle as pickle


def tex_f(f):
    return "${0:.3g}$".format(f)

data_path = "./data"
pb_table = load_pb_table(data_path)

f = open("go_id2name.p", 'rb')
go_id2name = pickle.load(f)
f.close()

f = open("go_id2cds_list.p", 'rb')
go_id2cds_list = pickle.load(f)
f.close()

table_tex_str = data_path + "/figure_2.tex"
table_tex = open(table_tex_str, 'w')
table_tex.writelines("\\documentclass[USLetter,5pt]{article}\n\\usepackage{adjustbox}\n")
table_tex.writelines("\\newcommand{\\specialcell}[2][c]{%\n\\begin{tabular}[#1]{@{}c@{}}#2\\end{tabular}}\n")
table_tex.writelines("\\begin{document}\n")

param = "siteomega-predmutsel"
label = "$\\omega_A =  \\omega - \\omega_0 $"
cds_id2omega = {}
table = str_to_table(pb_table, param)
for i, cds_id_str in enumerate(pb_table["CdsId"]):
    cds_id2omega[cds_id_str[2:-1]] = table[i]

p_value_list = []
mean_omega = np.mean(table)
std_omega = np.std(table)

outlier_set = set()
cds_file = open('{0}/outliers.out'.format(data_path), 'r')
cds_file.readline()
cds_file.readline()
for line in cds_file:
    if line != '\n':
        outlier_set.add(line.split("\t")[0])
cds_file.close()

cds_id_set = set(k for k, v in cds_id2omega.items())
cds_id_adap_set = outlier_set & cds_id_set
cds_id_no_adap_set = cds_id_set - cds_id_adap_set
print("Number of CDS under adaptation: {0}".format(len(cds_id_adap_set)))
for go_id, cds_list in go_id2cds_list.items():
    cds_go_set = set(cds_id for cds_id in cds_list)
    cds_no_go_set = cds_id_set - cds_go_set

    obs = np.array([[len(cds_go_set & cds_id_adap_set),
                     len(cds_no_go_set & cds_id_adap_set)],
                    [len(cds_go_set & cds_id_no_adap_set),
                     len(cds_no_go_set & cds_id_no_adap_set)]])
    assert sum(obs[0, :]) == len(cds_id_adap_set)
    assert sum(obs[1, :]) == len(cds_id_no_adap_set)
    assert sum(obs[:, 0]) == len(cds_go_set)
    assert sum(obs[:, 1]) == len(cds_no_go_set)
    assert np.sum(obs) == len(cds_id_set)

    if np.min(obs) > 1:
        oddsratio, p_value = st.fisher_exact(obs, alternative='greater')
        p_value_list.append((go_id, p_value, oddsratio, obs))

nbr_tests = len(p_value_list)
e_value_list = [p for p in p_value_list if nbr_tests * p[1] < 1]

heading = ["Gene Ontology", "$n_{\\mathrm{Observed}}$", "$n_{\\mathrm{Expected}}$", "Odds ratio", "$p_{\\mathrm{value}}$", "$e_{\\mathrm{value}}$"]

table_tex.writelines("\\begin{table}[ht]\n\\centering\n\\begin{adjustbox}{width = 1\\textwidth}\n\\small")
table_tex.writelines("\\begin{tabular}{" + "|c" * len(heading) + "|}\n")
table_tex.writelines("\\hline\n")
table_tex.writelines(" & ".join(heading) + "\\\\\n")
table_tex.writelines("\\hline\n")

e_value_list.sort(key=lambda x: x[1])

for go_id, p_value, oddsratio, obs in e_value_list:
    table_tex.writelines(" & ".join(
        [go_id2name[go_id].replace("_", "-"), tex_f(obs[0, 0]), tex_f(obs[0, 0] / oddsratio),
         tex_f(oddsratio), tex_f(p_value), tex_f(p_value * nbr_tests)]) + "\\\\\n")

caption = "\\textbf{" + label + "}. " + "{0} tests performed with {1} CDS detected as adaptive " \
                                        "and {2} as non-adaptive.\n".format(nbr_tests, len(cds_id_adap_set),
                                                                            len(cds_id_no_adap_set))
table_tex.writelines("\\hline\n")
table_tex.writelines("\\end{tabular}\n")
table_tex.writelines("\\end{adjustbox}\n" + "\\caption{" + caption + "}\n\\end{table}\n")

table_tex.writelines("\\end{document}\n")
table_tex.close()
print('Table generated')
os.system("pdflatex -output-directory={0} {1}".format(data_path, table_tex_str))
print('Pdf generated')

print('Ontology computed')
