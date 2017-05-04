import os
from lxml import etree

data_path = "/mnt/sda1/AdaptaPop/data"

homo_folder = "om_79_cds_homo"
homo_path = "{0}/{1}".format(data_path, homo_folder)

cds_file = open('{0}/cds_homo_pan.out'.format(data_path), 'w')
cds_file.write("Filename\tHomoCDS\tHomoTr\tPanCDS\tPanTr\n")

for file in os.listdir(homo_path):
    if file.endswith(".xml"):
        root = etree.parse("{0}/{1}.xml".format(homo_path, file[:-4])).getroot()
        cds_homo_79, tr_homo_79, cds_pan_79, tr_pan_79 = "", "", "", ""
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_homo_79 = specie.find('infoCDS').find('ensidTr').text
                cds_homo_79 = specie.find('infoCDS').find('ensid').text
            elif specie.attrib['species'] == 'Pan':
                tr_pan_79 = specie.find('infoCDS').find('ensidTr').text
                cds_pan_79 = specie.find('infoCDS').find('ensid').text
        cds_file.write("\t".join([file[:-4], cds_homo_79, tr_homo_79, cds_pan_79, tr_pan_79]) + "\n")

cds_file.close()

print("Done")