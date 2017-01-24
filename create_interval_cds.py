from lxml import etree
import os


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

hash_transcript = {}
gtf_file = open('Homo_sapiens_79_GRCh37.gtf', 'r')
hash_transcripts = {}
for line in gtf_file:
    line_split = line.split('\t')
    if len(line_split) > 7:
        info = line_split[8]
        if line_split[2] == 'CDS':
            transcript_find = info.find('transcript_id')
            exon = info[transcript_find + 15:transcript_find + 30]
            if not hash_transcripts.get(exon):
                hash_transcripts[exon] = Cds(line_split[0], line_split[6], exon)
            hash_transcripts[exon].add_exon(line_split[3], line_split[4])
gtf_file.close()

path = "om_79_cds_homo"
bedfile = open('79_interval_cds.bed', 'w')
for file in os.listdir("./" + path):
    if file.endswith(".xml"):
        root = etree.parse("./" + path + "/" + file).getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_id = specie.find('infoCDS').find('ensidTr').text
                for line in hash_transcripts[tr_id].befile_lines():
                    bedfile.write(line)
                break
bedfile.truncate()
bedfile.close()
