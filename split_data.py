import os
import numpy as np

path = "om_79_cds_mammals_no_pan_marsu"
tree = "mammals_no_homi_marsu.tree"
data = "/pandata/tlatrill/AdaptaPop/data"
# data = "./data"

alis = []
for file in os.listdir(data + "/" + path):
    if file.endswith("filtered_NT.ali"):
        alis.append(file)

assert len(alis) > 1, "Not enough files"
sample_1 = list(np.random.choice(alis, int(len(alis)/2), replace=False))
sample_2 = [ali for ali in alis if ali not in sample_1]

assert len(sample_1) + len(sample_2) == len(alis), "Sample not done correctly"

os.mkdir(data + "/" + path + "_1")
os.mkdir(data + "/" + path + "_2")
for s_1 in sample_1:
    os.system('cp ' + data + "/" + path + "/" + s_1 + " " + data + "/" + path + "_1")
for s_2 in sample_2:
    os.system('cp ' + data + "/" + path + "/" + s_2 + " " + data + "/" + path + "_2")

print('Spliting completed')
