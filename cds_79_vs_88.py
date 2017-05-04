import os
from Bio import SeqIO
from lxml import etree
from cds_libraries import build_dict_transcripts, build_dict_snps

data_path = "/mnt/sda1/AdaptaPop/data"
tmp_path = "/home/thibault/tmp"

homo_folder = "om_79_cds_homo"
homo_path = "{0}/{1}".format(data_path, homo_folder)
file_cds_homo_88 = "Homo_sapiens_88_GRCh38_cds_all"
file_cds_pan_88 = "Pan_88_cds_all"

tr_homo_79 = set()
cds_homo_79 = set()
tr_pan_79 = set()
cds_pan_79 = set()
for file in os.listdir(homo_path):
    if file.endswith(".xml"):
        root = etree.parse("{0}/{1}.xml".format(homo_path, file[:-4])).getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_homo_79.add(specie.find('infoCDS').find('ensidTr').text)
                cds_homo_79.add(specie.find('infoCDS').find('ensid').text)
            elif specie.attrib['species'] == 'Pan':
                tr_pan_79.add(specie.find('infoCDS').find('ensidTr').text)
                cds_pan_79.add(specie.find('infoCDS').find('ensid').text)
                break
print("Number of transcript from 79 Homo version: {0}".format(len(tr_homo_79)))
print("Number of cds from 79 Homo version: {0}".format(len(cds_homo_79)))
print("Number of transcript from 79 Pan version: {0}".format(len(tr_pan_79)))
print("Number of cds from 79 Pan version: {0}".format(len(cds_pan_79)))


tr_homo_88 = set()
cds_homo_88 = set()
for fasta in SeqIO.parse(open("{0}/{1}.fasta".format(data_path, file_cds_homo_88), 'r'), 'fasta'):
    tr_homo_88.add(fasta.id.split(".")[0])
    if "gene:" in fasta.description:
        cds_homo_88.add(fasta.description[fasta.description.find("gene")+len("gene:"):].split(" ")[0].split(".")[0])
print("Number of transcript from 88 Homo version: {0}".format(len(tr_homo_88)))
print("Number of cds from 88 Homo version: {0}".format(len(cds_homo_88)))

print("Number of transcript in the Homo intersection: {0}".format(len(tr_homo_79.intersection(tr_homo_88))))
print("Number of cds in the Homo intersection: {0}".format(len(cds_homo_79.intersection(cds_homo_88))))

tr_pan_88 = set()
cds_pan_88 = set()
for fasta in SeqIO.parse(open("{0}/{1}.fasta".format(data_path, file_cds_pan_88), 'r'), 'fasta'):
    tr_pan_88.add(fasta.id.split(".")[0])
    if "gene:" in fasta.description:
        cds_pan_88.add(fasta.description[fasta.description.find("gene")+len("gene:"):].split(" ")[0].split(".")[0])

print("Number of transcript from 88 Pan version: {0}".format(len(tr_pan_88)))
print("Number of cds from 88 Pan version: {0}".format(len(cds_pan_88)))
print("Number of transcript in the Pan intersection: {0}".format(len(tr_pan_79.intersection(tr_pan_88))))
print("Number of cds in the Pan intersection: {0}".format(len(cds_pan_79.intersection(cds_pan_88))))

print("Number of transcript in the all intersection: {0}".format(len(tr_homo_79.intersection(tr_homo_88).intersection(tr_pan_88))))
print("Number of cds in the all intersection: {0}".format(len(cds_homo_79.intersection(cds_homo_88).intersection(cds_pan_88))))

