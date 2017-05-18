import os
import numpy as np
from cds_libraries import load_snp_table

folders = ["siteomega", "mutsel", "mutselfreeomega"]

data_path = "/mnt/sda1/AdaptaPop/data"

cds_omega = {}
header = set()
header.add("pos")
for folder in folders:
    folder_path = "{0}/pb_cleaned_{1}".format(data_path, folder)
    for file in os.listdir(folder_path):
        if file.endswith(".predsiteomega") or file.endswith(".trace"):
            if file.endswith(".predsiteomega"):
                dtype = np.dtype([('#iter', 'str'), ('omega', 'float')])
                cols = [0, 1]
            elif folder == "siteomega":
                dtype = np.dtype([('#iter', 'int'), ('omega', 'float'), ('pos', 'float')])
                cols = [0, 6, 8]
            else:
                dtype = np.dtype([('#iter', 'int'), ('omega', 'float')])
                cols = [0, 6]
            table = np.loadtxt("{0}/{1}".format(folder_path, file), dtype=dtype, skiprows=1, usecols=cols)
            chain_name = file.replace('.predsiteomega', '').replace('.trace', '').replace('_filtered_NT', '')
            cds_name = chain_name[:chain_name.rfind("_")]
            if cds_name not in cds_omega:
                cds_omega[cds_name] = {}
            if file.endswith(".predsiteomega"):
                omega_name = "pred" + folder
            else:
                omega_name = folder
            if omega_name not in cds_omega[cds_name]:
                    cds_omega[cds_name][omega_name] = []
            header.add(omega_name)
            if file.endswith(".predsiteomega"):
                cds_omega[cds_name][omega_name].append(np.mean(table['omega']))
            else:
                cds_omega[cds_name][omega_name].append(np.mean(table['omega'][table['#iter'] > 100]))
                if folder == "siteomega":
                    if "pos" not in cds_omega[cds_name]:
                        cds_omega[cds_name]["pos"] = []
                    cds_omega[cds_name]["pos"].append(np.mean(table['pos'][table['#iter'] > 100]))

snp_table = load_snp_table(data_path, "88", "GRCh38")
snp_set = set([k for k, v in snp_table.items() if v["Chromosome"] not in ["X", "Y", "MT"]])

header = sorted(header)
txt_file = open(data_path + '/79_GRCh38_estimates_pb.out', 'w')
txt_file.write("CdsId\t"+"\t".join(header)+"\n")
for tr_id, omega_dico in cds_omega.items():
    line = tr_id
    if tr_id in snp_set:
        flag = True
    else:
        flag = False
    for omega_name in header:
        if omega_name in omega_dico:
            mean = np.mean(omega_dico[omega_name])
            if np.isnan(mean):
                flag = False
            line += "\t{0}".format(mean)
        else:
            flag = False
            line += "\tNaN"
    if flag:
        txt_file.write(line+"\n")
txt_file.close()
print("Job completed")
