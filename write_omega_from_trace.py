import os
import numpy as np

columns = ["globalomega", "siteomega", "mutsel", "mutselfreeomega"]

# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

cds_omega = {}
header = set()
for column in columns:
    folder_path = "{0}/pb_{1}".format(data_path, column)
    for file in os.listdir(folder_path):
        if file.endswith(".predsiteomega") or file.endswith(".trace"):
            if file.endswith(".predsiteomega"):
                dtype = np.dtype([('#iter', 'str'), ('omega', 'float')])
                cols = [0, 1]
            else:
                dtype = np.dtype([('#iter', 'int'), ('omega', 'float')])
                cols = [0, 6]
            table = np.loadtxt("{0}/{1}".format(folder_path, file), dtype=dtype, skiprows=1, usecols=cols)
            chain_name = file.replace('.predsiteomega', '').replace('.trace', '').replace('_filtered_NT', '')
            cds_name = chain_name[:chain_name.rfind("_")]
            if cds_name not in cds_omega:
                cds_omega[cds_name] = {}
            if file.endswith(".predsiteomega"):
                omega_name = "pred" + column
            else:
                omega_name = column
            if omega_name not in cds_omega[cds_name]:
                    cds_omega[cds_name][omega_name] = []
            header.add(omega_name)
            if file.endswith(".predsiteomega"):
                cds_omega[cds_name][omega_name].append(np.mean(table['omega']))
            else:
                cds_omega[cds_name][omega_name].append(np.mean(table['omega'][table['#iter'] > 100]))

header = sorted(header)
txt_file = open(data_path + '/79_omega_estimated.out', 'w')
txt_file.write("CdsId\t"+"\t".join(header)+"\n")
for tr_id, omega_dico in cds_omega.items():
    line = tr_id
    for omega_name in header:
        if omega_name in omega_dico:
            line += "\t{0}".format(np.mean(omega_dico[omega_name]))
        else:
            line += "\tnan"
    txt_file.write(line+"\n")
txt_file.close()
print("Job completed")
