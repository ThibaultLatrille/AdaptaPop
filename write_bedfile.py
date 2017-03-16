from lxml import etree
import os

# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"


class Cds(object):
    def __init__(self, chromosome, strand, name):
        self.chromosome = chromosome
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        assert int(start_exon) <= int(end_exon), 'The start position is greater than the end position of the exon '
        self.exons.append((int(start_exon), int(end_exon)))
        if len(self.exons) > 1 and self.exons[-1][0] < self.exons[-2][1]:
            self.exons.sort(key=lambda x: x[0])
            for i in range(len(self.exons) - 1):
                assert self.exons[i][1] < self.exons[i + 1][0], "At least one exon is overlapping with an other"

    def befile_lines(self):
        lines = []
        for start, end in self.exons:
            lines.append(self.chromosome + "\t" + str(start) + "\t" + str(end) +
                         "\t" + self.name + "\t0\t" + self.strand + "\n")
        return lines


def build_dict_transcripts(in_data, file_name):
    gtf_file = open(in_data + "/" + file_name, 'r')
    in_dict_transcripts = {}
    in_not_confirmed_cds = {}
    for line in gtf_file:
        line_split = line.split('\t')
        if len(line_split) > 7:
            info = line_split[8]
            if line_split[2] == 'CDS':
                transcript_find = info.find('transcript_id')
                if transcript_find != -1:
                    tr_id = info[transcript_find + 15:transcript_find + 30]
                    if info.find('cds_start_NF') != -1 or info.find('cds_end_NF') != -1:
                        if not in_not_confirmed_cds.get(tr_id):
                            in_not_confirmed_cds[tr_id] = True
                    if not in_dict_transcripts.get(tr_id):
                        in_dict_transcripts[tr_id] = Cds(line_split[0], line_split[6], tr_id)
                        in_dict_transcripts[tr_id].add_exon(line_split[3], line_split[4])
    gtf_file.close()
    return in_dict_transcripts, in_not_confirmed_cds


def build_bedfile(in_data, path, file_name, in_dict_transcripts):
    in_bedfile = open(in_data + "/" + file_name, 'w')
    for file in os.listdir(in_data + "/" + path):
        if file.endswith(".xml"):
            root = etree.parse(in_data + "/" + path + "/" + file).getroot()
            for specie in root.find('CDS').findall("speciesCDS"):
                if specie.attrib['species'] == 'Homo':
                    tr_id = specie.find('infoCDS').find('ensidTr').text
                    if in_dict_transcripts.get(tr_id):
                        for line in in_dict_transcripts[tr_id].befile_lines():
                            in_bedfile.write(line)
                        break
    in_bedfile.truncate()
    in_bedfile.close()


dict_transcripts, not_confirmed_cds = build_dict_transcripts(data, 'Homo_sapiens_79_GRCh37.gtf')
build_bedfile(data, "om_79_cds_homo", '79_interval_cds.bed', dict_transcripts)

print('Job completed')