import os
from Bio import SeqIO
from lxml import etree
from collections import defaultdict

# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"


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

    def nt_position(self, position):
        nt = -1
        running_sum = 0
        if self.strand == "+":
            for start, end in self.exons:
                if start <= position <= end:
                    nt = position - start + running_sum
                    break
                running_sum += end - start + 1
        else:
            for start, end in reversed(self.exons):
                if start <= position <= end:
                    nt = end - position + running_sum
                    break
                running_sum += end - start + 1
        assert nt > -1, "The position given is not inside the cds"
        return nt

    def amino_acid(self, seq, position, nt_ref, nt_alt):
        nt_position = self.nt_position(position)
        frame = nt_position % 3
        codon_ref = str(seq[nt_position - frame:nt_position + 3 - frame])
        aa_ref = codontable[codon_ref]
        if self.strand == "-":
            nt_ref = complement[nt_ref]
            nt_alt = complement[nt_alt]
        codon_alt = codon_ref[:frame] + nt_alt + codon_ref[frame + 1:]
        aa_alt = codontable[codon_alt]
        if codon_ref[frame] != nt_ref:
            return '!', codon_ref[frame], codon_ref, codon_alt
        else:
            return aa_ref, aa_alt, codon_ref, codon_alt

    def seq_length(self):
        return sum([j - i + 1 for i, j in self.exons])


def build_hash_transcripts(data, file_name):
    gtf_file = open(data + "/" + file_name, 'r')
    hash_transcripts = {}
    not_confirmed_cds = {}
    for line in gtf_file:
        line_split = line.split('\t')
        if len(line_split) > 7:
            info = line_split[8]
            if line_split[2] == 'CDS':
                transcript_find = info.find('transcript_id')
                if transcript_find != -1:
                    tr_id = info[transcript_find + 15:transcript_find + 30]
                    if info.find('cds_start_NF') != -1 or info.find('cds_end_NF') != -1:
                        if not not_confirmed_cds.get(tr_id):
                            not_confirmed_cds[tr_id] = True
                    if not hash_transcripts.get(tr_id):
                        hash_transcripts[tr_id] = Cds(line_split[0], line_split[6], tr_id)
                    hash_transcripts[tr_id].add_exon(line_split[3], line_split[4])
    gtf_file.close()
    return hash_transcripts, not_confirmed_cds


def build_hash_snps(data, file_name):
    vcf_file = open(data + "/" + file_name, 'r')
    hash_snps = {}
    for line in vcf_file:
        if line[0] != '#':
            split_line = line.split("\t")
            if len(split_line[3]) == 1 and len(split_line[4]) == 1 and split_line[4] != ".":
                tr_id = split_line[11]
                if not hash_snps.get(tr_id):
                    hash_snps[tr_id] = []
                hash_snps[tr_id].append((split_line[2], split_line[0], int(split_line[1]), split_line[3], split_line[4]))
    vcf_file.close()
    return hash_snps


homo_seq_ungap, homo_seq, pan_seq, tr_id = "", "", "", ""
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'stop', 'TAG': 'stop',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'stop', 'TGG': 'W'})

hash_transcripts, not_confirmed_cds = build_hash_transcripts(data, 'Homo_sapiens_79_GRCh37.gtf')
hash_snps = build_hash_snps(data, 'Homo_sapiens_79_polymorphism_in_cds.vcf')

path = "om_79_cds_homo"
txt_file = open('79_mk_test.out', 'w')
txt_file.write("TrId\tSeqLength\tNbrExons\tPtot\tPn\tPs\tDtot\tDn\tDs\n")
cds_total = 0
errors_cds_nf = []
errors_cds_length = []
errors_cds_ungap_length = []
errors_cds_unequal_length = []
snp_total = 0
errors_snp_ref = []
errors_snp_stop = []
errors_snp_alt_stop = []
errors_snp_codon = []
for file in os.listdir(data + "/" + path):
    if file.endswith(".xml"):
        cds_total += 1
        file_name = file[:-4]

        root = etree.parse(data + "/" + path + "/" + file_name + ".xml").getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_id = specie.find('infoCDS').find('ensidTr').text
                break

        if not_confirmed_cds.get(tr_id):
            errors_cds_nf.append(tr_id)
            continue

        cds = hash_transcripts[tr_id]
        cds_length = cds.seq_length()
        if cds_length % 3 != 0:
            errors_cds_length.append(tr_id)
            continue

        fasta_seqs = SeqIO.parse(open(data + "/" + path + "/" + file_name + "_raw_NT.fasta", 'r'), 'fasta')
        for fasta in fasta_seqs:
            if fasta.id == "Homo":
                homo_seq = fasta.seq
                homo_seq_ungap = homo_seq.ungap('-')
            elif fasta.id == "Pan":
                pan_seq = fasta.seq

        ungap_length = len(homo_seq_ungap)
        if ungap_length % 3 != 0:
            errors_cds_ungap_length.append(tr_id)
            continue

        if cds_length != ungap_length - 3:
            errors_cds_unequal_length.append(tr_id)
            continue

        if hash_snps.get(tr_id):
            list_aa_poly = []
            for snp_id, chromosome, pos, ref_nt, alt_nt in hash_snps[tr_id]:
                snp_total += 1
                ref_aa, alt_aa, ref_codon, alt_codon = cds.amino_acid(homo_seq_ungap, pos, ref_nt, alt_nt)
                if ref_aa == '!':
                    errors_snp_ref.append((tr_id, snp_id, chromosome, str(pos), alt_aa, ref_nt, alt_nt, ref_codon, alt_codon))
                elif ref_aa == '-' or alt_aa == '-':
                    errors_snp_codon.append((tr_id, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif ref_aa == 'stop':
                    errors_snp_stop.append((tr_id, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif alt_aa == 'stop':
                    errors_snp_alt_stop.append((tr_id, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                else:
                    list_aa_poly.append((ref_aa, alt_aa))

        list_aa_div = []
        for pos, (homo_n, pan_n) in enumerate(zip(homo_seq, pan_seq)):
            if homo_n != pan_n and homo_n != '-' and pan_n != '-':
                homo_nt_pos = len(homo_seq[:pos].ungap('-'))
                homo_frame = homo_nt_pos % 3
                homo_codon = str(homo_seq_ungap[homo_nt_pos - homo_frame:homo_nt_pos + 3 - homo_frame])
                pan_nt_pos = len(pan_seq[:pos].ungap('-'))
                pan_frame = pan_nt_pos % 3
                pan_codon = str(pan_seq.ungap('-')[pan_nt_pos - pan_frame:pan_nt_pos + 3 - pan_frame])
                list_aa_div.append((codontable[homo_codon], codontable[pan_codon]))

        pn = str(len([0 for i, j in list_aa_poly if i != j]))
        ps = str(len([0 for i, j in list_aa_poly if i == j]))
        dn = str(len([0 for i, j in list_aa_div if i != j and i != '-' and j != '-']))
        ds = str(len([0 for i, j in list_aa_div if i == j and i != '-' and j != '-']))
        txt_file.write(file_name + "\t" + str(len(homo_seq_ungap)) + "\t" + str(len(cds.exons)) + "\t"
                       + str(len(list_aa_poly)) + "\t" + pn + "\t" + ps + "\t"
                       + str(len(list_aa_div)) + "\t" + dn + "\t" + ds + "\n")
txt_file.close()

error_file = open('79_mk_test_errors.out', 'w')
cds_errors = len(errors_cds_nf)+len(errors_cds_length)+len(errors_cds_ungap_length)+len(errors_cds_unequal_length)
if cds_errors > 0:
    error_file.write(str(cds_errors) + " errors out of " + str(cds_total) + " coding sequences ("
                     + '%.3f' % (cds_errors*100./cds_total) + "%)")
    if len(errors_cds_nf) > 0:
        error_file.write("\n\nCoding region start or end could not be confirmed ("
                         + str(len(errors_cds_nf)) + " cds):\n")
        error_file.write(" ".join(errors_cds_nf))
    if len(errors_cds_length) > 0:
        error_file.write("\n\nThe computed sequence size from exons starts and ends is not a multiple of 3 ("
                         + str(len(errors_cds_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_length))
    if len(errors_cds_ungap_length) > 0:
        error_file.write("\n\nThe sequence size is not a multiple of 3 ("
                         + str(len(errors_cds_ungap_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_ungap_length))
    if len(errors_cds_unequal_length) > 0:
        error_file.write("\n\nThe computed sequence size from exons starts and ends doesn't match the sequence size ("
                         + str(len(errors_cds_unequal_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_unequal_length))

snp_errors = len(errors_snp_ref)+len(errors_snp_codon)+len(errors_snp_stop)+len(errors_snp_alt_stop)
if snp_errors > 0:
    error_file.write("\n\n" + str(snp_errors) + " errors out of " + str(snp_total) + " SNPs ("
                     + '%.3f' % (snp_errors*100./snp_total) + "%)")
    if len(errors_snp_ref) > 0:
        error_file.write("\n\nSNPs retrieved from the fasta are not equal to the reference ("
                         + str(len(errors_snp_ref)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tSeqNt\tRefNt\tAltNt\tRefCodon\tAltCodon\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_ref]))
    if len(errors_snp_codon) > 0:
        error_file.write("\n\n The reference or alternate amino-acid can't be identified ("
                         + str(len(errors_snp_codon)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_codon]))
    if len(errors_snp_stop) > 0:
        error_file.write("\n\n The reference amino-acid is a stop codon ("
                         + str(len(errors_snp_stop)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_stop]))
    if len(errors_snp_alt_stop) > 0:
        error_file.write("\n\n The alternate amino-acid is a stop codon ("
                         + str(len(errors_snp_alt_stop)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_alt_stop]))
error_file.close()

print("Job completed")
