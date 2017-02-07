from lxml import etree
import os
import matplotlib.pyplot as plt
import numpy as np

mk_file = "79_mk_test.out"
pbmpi_path = "pb_output"
# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"

mk_dict = {}
pb_dict = {}
nbr_bins = 5

mk_data = open(data + "/" + mk_file, 'r')
mk_header = mk_data.readline().split('\t')
for line in mk_data:
    line_split = line.split('\t')
    mk_dict[line_split[0]] = dict(zip(mk_header[1:], [float(i) for i in line_split[1:]]))
mk_data.close()


def alpha(in_data, path, specie_name):
    dict_cds_to_tr = {}
    for in_file in os.listdir(in_data + "/" + path):
        if in_file.endswith(".xml"):
            root = etree.parse(data + "/" + path + "/" + in_file).getroot()
            for specie in root.find('CDS').findall("speciesCDS"):
                if specie.attrib['species'] == specie_name:
                    tr_id = specie.find('infoCDS').find('ensidTr').text
                    dict_cds_to_tr[in_file[:-4]] = tr_id
    return dict_cds_to_tr

for file in os.listdir(data + "/" + pbmpi_path):
    if file.endswith(".trace"):
        pb_data = open(data + "/" + pbmpi_path + "/" + file, 'r')
        chain_name = file.replace('.trace', '').replace('_filtered_NT', '')
        cds_name = chain_name[:chain_name.rfind("_")]
        line = pb_data.readline()
        pb_header = line.split('\t')
        for line in pb_data:
            pass
        if pb_dict.get(cds_name):
            for key, value in zip(pb_header, line.split("\t")):
                pb_dict[cds_name][key].append(float(value))
        else:
            pb_dict[cds_name] = dict(zip(pb_header, [[float(i)] for i in line.split("\t")]))
        pb_data.close()

pb_mean_omega = [(tr_id, np.mean(trace_dict["omega"])) for tr_id, trace_dict in pb_dict.items()]
pb_mean_omega.sort(key=lambda x: x[1])

bin_size = int(int(len(pb_mean_omega)) / nbr_bins)
grouped_omega = [pb_mean_omega[i:i + bin_size] for i in range(0, len(pb_mean_omega), bin_size)]

grouped_alpha = [[mk_dict[tr_id]["Pn"] for tr_id, _ in group] for group in grouped_omega]
grouped_omega = [[omega for _, omega in group] for group in grouped_omega]

list_alpha = [np.mean(group) for group in grouped_alpha]
list_omega = [np.mean(group) for group in grouped_omega]

plt.plot(list_omega, list_alpha, linewidth=3)
print(list_alpha)
print(list_omega)
plt.show()
print('Test completed')
