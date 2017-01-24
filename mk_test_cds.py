import os
from Bio import SeqIO
from lxml import etree
from collections import defaultdict


class Cds(object):
    def __init__(self, chromosome, strand, name):
        self.chromosome = chromosome
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        assert int(start_exon) <= int(end_exon), 'The start position is greater than the end position of the exon'
        self.exons.append((int(start_exon), int(end_exon)))
        if len(self.exons) > 1 and self.exons[-1][0] < self.exons[-2][1]:
            self.exons.sort(key=lambda x: x[0])
            for i in range(len(self.exons) - 1):
                assert self.exons[i][1] < self.exons[i + 1][0], "At least one exon is overlapping with an other"

    def nucleotide_position(self, position):
        nucleotide = -1
        running_sum = 0
        if self.strand == "+":
            for start, end in self.exons:
                if start <= position <= end:
                    nucleotide = position - start + running_sum
                    break
                running_sum += end - start + 1
        else:
            for start, end in reversed(self.exons):
                if start <= position <= end:
                    nucleotide = end - position + running_sum
                    break
                running_sum += end - start + 1
        assert nucleotide > -1, "The position given is not inside the cds"
        return nucleotide

    def amino_acid(self, seq, position, ref_nucleotide, alt_nucleotide):
        nucleotide_position = self.nucleotide_position(position)
        frame = nucleotide_position % 3
        codon = str(seq[nucleotide_position - frame:nucleotide_position + 3 - frame])
        ref_aa = codontable[codon]
        if self.strand == "-":
            ref_nucleotide = complement[ref_nucleotide]
            alt_nucleotide = complement[alt_nucleotide]
        # assert codon[frame] == ref_nucleotide, 'Reference nucleotide and retrieved nucleotide ' \
        #                                        'from the reference sequence are not the same'
        alt_aa = codontable[codon[:frame] + alt_nucleotide + codon[frame + 1:]]
        return ref_aa, alt_aa

    def sequence_length(self):
        return sum([j - i + 1 for i, j in self.exons])


start, end = 0, 0
homo_sequence_ungap, homo_sequence, pan_sequence, tr_id = "", "", "", ""
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
codontable = defaultdict(lambda: "-")
codontable.update()

os.chdir('/home/thibault/adaptapop')
gtf_file = open('Homo_sapiens_79_GRCh37.gtf', 'r')
hash_transcripts = {}
for line in gtf_file:
    line_split = line.split('\t')
    if len(line_split) > 7:
        info = line_split[8]
        if line_split[2] == 'CDS':
            transcript_find = info.find('transcript_id')
            tr_id = info[transcript_find + 15:transcript_find + 30]
            if not hash_transcripts.get(tr_id):
                hash_transcripts[tr_id] = Cds(line_split[0], line_split[6], tr_id)
            hash_transcripts[tr_id].add_exon(line_split[3], line_split[4])
gtf_file.close()

os.chdir('/home/thibault/adaptapop')
vcf_file = open('Homo_sapiens_79_polymorphism_in_cds.vcf', 'r')
hash_variations = {}
for line in vcf_file:
    if line[0] != '#':
        split_line = line.split("\t")
        if len(split_line[3]) == 1 and len(split_line[4]) == 1 and split_line[4] != ".":
            tr_id = split_line[11]
            if not hash_variations.get(tr_id):
                hash_variations[tr_id] = []
            hash_variations[tr_id].append((int(split_line[1]), split_line[3], split_line[4]))
vcf_file.close()

path = "om_79_cds_homo"
txt_file = open('79_mk_test.out', 'w')
for file in os.listdir("./" + path):
    if file.endswith(".xml"):
        file_name = file[:-4]
        txt_file.write("\nTranscript: " + file_name)
        root = etree.parse("./" + path + "/" + file_name + ".xml").getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_id = specie.find('infoCDS').find('ensidTr').text
                break

        fasta_sequences = SeqIO.parse(open("./" + path + "/" + file_name + "_raw_NT.fasta", 'r'), 'fasta')
        for fasta in fasta_sequences:
            if fasta.id == "Homo":
                homo_sequence = fasta.seq
                homo_sequence_ungap = homo_sequence.ungap('-')
            elif fasta.id == "Pan":
                pan_sequence = fasta.seq

        cds = hash_transcripts[tr_id]
        txt_file.write("\nSequence length:" + str(len(homo_sequence_ungap)))
        txt_file.write("\nNumber of exons:" + str(len(cds.exons)))

        if hash_variations.get(tr_id):
            list_aa_poly = []
            for position, ref_nucleotide, alt_nucleotide in hash_variations[tr_id]:
                list_aa_poly.append(cds.amino_acid(homo_sequence_ungap, position, ref_nucleotide, alt_nucleotide))

            txt_file.write("\nNumber of polymorphic SNPs:" + str(len(list_aa_poly)))
            txt_file.write("\nNumber of polymorphic synonymous SNPs:" + str(len([0 for i, j in list_aa_poly
                                                                                 if i == j and i != '-' and j != '-'])))
            txt_file.write("\nNumber of polymorphic non-synonymous SNPs:" + str(len([0 for i, j in list_aa_poly
                                                                                     if i != j and i != '-' and j != '-'])))

        list_div = []
        for position, (homo_n, pan_n) in enumerate(zip(homo_sequence, pan_sequence)):
            if homo_n != pan_n and homo_n != '-' and pan_n != '-':
                list_div.append((position, homo_n, pan_n))

        list_aa_div = []
        for position, homo_n, pan_n in list_div:
            homo_nucleotide_position = len(homo_sequence[:position].ungap('-'))
            homo_frame = homo_nucleotide_position % 3
            homo_codon = str(
                homo_sequence_ungap[homo_nucleotide_position - homo_frame:homo_nucleotide_position + 3 - homo_frame])
            pan_nucleotide_position = len(pan_sequence[:position].ungap('-'))
            pan_frame = pan_nucleotide_position % 3
            pan_codon = str(
                pan_sequence.ungap('-')[pan_nucleotide_position - pan_frame:pan_nucleotide_position + 3 - pan_frame])
            list_aa_div.append((codontable[homo_codon], codontable[pan_codon]))

        txt_file.write("\nNumber of divergence:" + str(len(list_div)))
        txt_file.write("\nNumber of synonymous divergence:" + str(len([0 for i, j in list_aa_div
                                                                       if i == j and i != '-' and j != '-'])))
        txt_file.write("\nNumber of non-synonymous divergence:" + str(len([0 for i, j in list_aa_div
                                                                           if i != j and i != '-' and j != '-'])))
        txt_file.write("\n")
txt_file.close()
