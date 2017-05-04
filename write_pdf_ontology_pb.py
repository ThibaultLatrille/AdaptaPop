import numpy as np
import os
from lxml import etree
import scipy.stats as st
from cds_libraries import params_pb, load_pb_table, str_to_table


def tex_f(f):
    return "${0:.3g}$".format(f)


cds_folder = "om_79_cds_homo"
data_path = "/mnt/sda1/AdaptaPop/data"
cds_path = "{0}/{1}".format(data_path, cds_folder)
pb_table = load_pb_table(data_path)

pb_cds_ids = set(cds_id[2:-1] for cds_id in pb_table["CdsId"])
go_id2name = {}
go_id2cds_list = {}

cds_list = open('{0}/cds_homo_pan.out'.format(data_path), 'r')
cds_list_header = cds_list.readline().replace('\n', '').split('\t')
for line in cds_list:
    line_split = line.replace('\n', '').split('\t')
    cds_dict = dict(zip(cds_list_header, line_split))
    cds_id = cds_dict["HomoCDS"]
    file_name = cds_dict["Filename"]
    if cds_id in pb_cds_ids:
        root = etree.parse("{0}/{1}.xml".format(cds_path, file_name)).getroot()
        for annot in root.find('goAnnot').findall("annot"):
            go_id = annot.find('goId').text
            go_name = annot.find('goName').text
            if go_id not in go_id2name:
                go_id2name[go_id] = go_name
            if go_id not in go_id2cds_list:
                go_id2cds_list[go_id] = [cds_id]
            else:
                go_id2cds_list[go_id].append(cds_id)

parametric = False
for negative in [True, False]:
    table_tex_str = data_path + "/79_GRCh38_ontology_table_{0}{1}.tex".format("parametric" if parametric else "wilcoxon", "" if negative else "_only_pos")
    table_tex = open(table_tex_str, 'w')
    table_tex.writelines("\\documentclass[USLetter,5pt]{article}\n\\usepackage{adjustbox}\n")
    table_tex.writelines("\\newcommand{\\specialcell}[2][c]{%\n\\begin{tabular}[#1]{@{}c@{}}#2\\end{tabular}}\n")
    table_tex.writelines("\\begin{document}\n")

    txt_file = open(data_path + '/79_GRCh38_ontology_{0}{1}.out'.format("parametric" if parametric else "wilcoxon", "" if negative else "_only_pos"), 'w')
    for param, label in params_pb:
        cds_id2omega = {}
        table = str_to_table(pb_table, param)
        for i, cds_id_str in enumerate(pb_table["CdsId"]):
            cds_id2omega[cds_id_str[2:-1]] = table[i]

        p_value_list = []
        table = [val for val in table if negative or val > 0]
        mean_omega = np.mean(table)
        std_omega = np.std(table)

        cds_id_set = set([k for k, v in cds_id2omega.items() if negative or v > 0])
        for go_id, cds_list in go_id2cds_list.items():
            cds_list_filter = [cds_id for cds_id in cds_list if negative or cds_id2omega[cds_id] > 0]
            n = len(cds_list_filter)
            if n > 20:
                if parametric:
                    go_mean_omega = np.mean([cds_id2omega[cds_id] for cds_id in cds_list_filter])
                    score = np.sqrt(n) * (go_mean_omega - mean_omega) / std_omega
                    p_value = 1 - st.norm.cdf(score)
                    p_value_list.append((go_id, p_value, score, n, go_mean_omega))
                else:
                    list_in = [cds_id2omega[cds_id] for cds_id in cds_list_filter]
                    list_out = [cds_id2omega[cds_id] for cds_id in cds_id_set.difference(set(cds_list_filter))]
                    assert len(list_in) + len(list_out) == len(table)
                    score, p_value = st.mannwhitneyu(list_in, list_out, alternative="greater")
                    p_value_list.append((go_id, p_value, score, n, np.mean(list_in)))

        nbr_cds = len(p_value_list)
        e_value_list = [p for p in p_value_list if nbr_cds * p[1] < 1]

        if len(e_value_list) > 0:
            heading = ["goId", "goName", "n", "p-value", "e-value", "z-score" if parametric else "Mann-Whitney-U", label]

            txt_file.write("nbr of CDS evaluated = {0}\n".format(nbr_cds))
            txt_file.write("mean {0} = {1:.3g}\n".format(label, mean_omega))
            txt_file.write("std {0} = {1:.3g}\n".format(label, std_omega))
            txt_file.write("\t".join(heading) + "\n")

            table_tex.writelines("\\begin{table}[ht]\n\\centering\n\\begin{adjustbox}{width = 1\\textwidth}\n\\small")
            table_tex.writelines("\\begin{tabular}{" + "|c" * len(heading) + "|}\n")
            table_tex.writelines("\\hline\n")
            table_tex.writelines(" & ".join(heading) + "\\\\\n")
            table_tex.writelines("\\hline\n")

            e_value_list.sort(key=lambda x: x[1])

            for go_id, p_value, score, n, go_mean_omega in e_value_list:
                txt_file.write(
                    "{0}\t{1}\t{2}\t{3:.3g}\t{4:.3g}\t{5:.3g}\t{6:.3g}\n".format(go_id, go_id2name[go_id], n, p_value,
                                                                                 p_value * nbr_cds, score,
                                                                                 go_mean_omega))
            txt_file.write("\n")

            for go_id, p_value, score, n, go_mean_omega in e_value_list:
                table_tex.writelines(" & ".join(
                    [go_id, go_id2name[go_id].replace("_", "-"), tex_f(n), tex_f(p_value), tex_f(p_value * nbr_cds), tex_f(score),
                     tex_f(go_mean_omega)]) + "\\\\\n")

            caption = "\\textbf{" + label + "}. " + "{0} tests performed,  \n".format(nbr_cds)
            caption += "$E[${0}$]={1:.3g}$ and \n".format(label, mean_omega)
            caption += "$\\sigma(${0}$)={1:.3g}$\n".format(label, std_omega)
            table_tex.writelines("\\hline\n")
            table_tex.writelines("\\end{tabular}\n")
            table_tex.writelines("\\end{adjustbox}\n" + "\\caption{" + caption + "}\n\\end{table}\n")

    table_tex.writelines("\\end{document}\n")
    txt_file.close()
    table_tex.close()
    print('Table generated')
    os.system("pdflatex -output-directory={0} {1}".format(data_path, table_tex_str))
    print('Pdf generated')

print('Ontology computed')
