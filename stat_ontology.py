import os
import matplotlib as mpl
import numpy as np
from lxml import etree
import scipy.stats as st

cds_folder = "om_79_cds_mammals_no_pan_marsu"
omega_estimated = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"
cds_path = "{0}/{1}".format(data_path, cds_folder)

dtype = np.dtype([("tr_id", 'str', 32)] + [(name, 'float64', 1) for name in omega_estimated])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)

params = [("siteomega", "omega"),
          ("predmutsel", "omega_0"),
          ("predmutselfreeomega", "omega_0^*"),
          ("mutselfreeomega", "omega_*"),
          ("siteomega/predmutsel", "omega/omega_0"),
          ("siteomega-predmutsel", "omega_a"),
          ("predmutselfreeomega*(mutselfreeomega-1)", "omega_a^*")]


def str_to_table(table, label, f_list):
    return eval(label, {f: table[f] for f in f_list})

go_id2name = {}
go_id2cds_list = {}
for cds_id_str in omega_table["tr_id"]:
    cds_id = cds_id_str[2:-1]
    root = etree.parse("{0}/{1}.xml".format(cds_path, cds_id)).getroot()
    for annot in root.find('goAnnot').findall("annot"):
        go_id = annot.find('goId').text
        go_name = annot.find('goName').text
        if go_id not in go_id2name:
            go_id2name[go_id] = go_name
        if go_id not in go_id2cds_list:
            go_id2cds_list[go_id] = [cds_id]
        else:
            go_id2cds_list[go_id].append(cds_id)


txt_file = open(data_path + '/79_ontology.out', 'w')
for omega, name in params:
    cds_id2omega = {}
    table = str_to_table(omega_table, omega, omega_estimated)
    for i, cds_id_str in enumerate(omega_table["tr_id"]):
        cds_id2omega[cds_id_str[2:-1]] = table[i]

    p_value_list = []

    mean_omega = np.mean(table)
    std_omega = np.std(table)
    for go_id, cds_list in go_id2cds_list.items():
        n = len(cds_list)
        if n > 4:
            go_mean_omega = np.mean([cds_id2omega[cds_id] for cds_id in cds_list])
            z_score = np.sqrt(n)*(go_mean_omega - mean_omega) / std_omega
            p_value = 1 - st.norm.cdf(z_score)
            p_value_list.append((go_id, p_value, z_score, n, go_mean_omega))

    nbr_ontologies = len(p_value_list)
    p_value_list = [p for p in p_value_list if p[1] < 0.05 / nbr_ontologies]

    if len(p_value_list) > 0:
        txt_file.write("mean {0} = {1:.3g}\n".format(name, mean_omega))
        txt_file.write("std {0} = {1:.3g}\n".format(name, std_omega))
        txt_file.write("goId\tgoName\tn\tp-value\tz-score\t{0}\n".format(name))
        p_value_list.sort(key=lambda x: x[1])
        for go_id, p_value, z_score, n, go_mean_omega in p_value_list:
            txt_file.write("{0}\t{1}\t{2}\t{3:.3g}\t{4:.3g}\t{5:.3g}\n".format(go_id, go_id2name[go_id], n, p_value, z_score,  go_mean_omega))
        txt_file.write("\n")
txt_file.close()
print('Ontology computed')
