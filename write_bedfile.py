from lxml import etree
import os
from cds_libraries import build_dict_transcripts

# data = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"


def build_bedfile(_data_path, cds_folder, file_name, in_dict_transcripts):
    in_bedfile = open("{0}/{1}".format(_data_path,  file_name), 'w')
    for file in os.listdir("{0}/{1}".format(_data_path,  cds_folder)):
        if file.endswith(".xml"):
            root = etree.parse("{0}/{1}/[2}".format(_data_path,  cds_folder, file)).getroot()
            for specie in root.find('CDS').findall("speciesCDS"):
                if specie.attrib['species'] == 'Homo':
                    tr_id = specie.find('infoCDS').find('ensidTr').text
                    if in_dict_transcripts.get(tr_id):
                        for line in in_dict_transcripts[tr_id].befile_lines():
                            in_bedfile.write(line)
                        break
    in_bedfile.truncate()
    in_bedfile.close()


dict_transcripts, _, _ = build_dict_transcripts(data_path, 'Homo_sapiens_79_GRCh37.gtf')
build_bedfile(data_path, "om_79_cds_homo", '79_interval_cds.bed', dict_transcripts)

print('Job completed')
