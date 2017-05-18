from lxml import etree
from cds_libraries import load_pb_table
import pickle as pickle


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

file = open("go_id2cds_list.p", 'wb')
pickle.dump(go_id2cds_list, file)
file.close()

file = open("go_id2name.p", 'wb')
pickle.dump(go_id2name, file)
file.close()