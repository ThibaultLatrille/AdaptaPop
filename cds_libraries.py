from collections import defaultdict
import numpy as np


class Cds(object):
    def __init__(self, chromosome, strand, name):
        self.chromosome = chromosome
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        assert int(start_exon) <= int(end_exon), 'The start position is greater than the end position of the exon'
        if self.strand == "+":
            self.exons.append((int(start_exon), int(end_exon)))
        else:
            self.exons.insert(0, (int(start_exon), int(end_exon)))
        for i in range(len(self.exons) - 1):
            if not self.exons[i][1] < self.exons[i + 1][0]:
                print("At least one exon is overlapping with an other")

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

    def befile_lines(self):
        lines = []
        for start, end in self.exons:
            lines.append("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(self.chromosome, start, end, self.name, self.strand))
        return lines


def build_dict_transcripts(_data_path, file_name):
    gtf_file = open("{0}/{1}".format(_data_path, file_name), 'r')
    _dict_transcripts = {}
    _not_confirmed_cds = {}
    _full_transcripts = {}
    for line in gtf_file:
        line_split = line.split('\t')
        if len(line_split) > 7:
            info = line_split[8]
            if line_split[2] == 'CDS':
                transcript_find = info.find('transcript_id')
                if transcript_find != -1:
                    tr_id = info[transcript_find + 15:transcript_find + 30]
                    if info.find('cds_start_NF') != -1 or info.find('cds_end_NF') != -1:
                        if tr_id not in _not_confirmed_cds:
                            _not_confirmed_cds[tr_id] = True
                    if tr_id not in _dict_transcripts:
                        _dict_transcripts[tr_id] = Cds(line_split[0], line_split[6], tr_id)
                        _full_transcripts[tr_id] = [line]
                    else:
                        _full_transcripts[tr_id].append(line)
                    _dict_transcripts[tr_id].add_exon(line_split[3], line_split[4])

    gtf_file.close()
    return _dict_transcripts, _not_confirmed_cds, _full_transcripts


def build_dict_snps(_data_path, file_name):
    vcf_file = open("{0}/{1}".format(_data_path, file_name), 'r')
    _dict_snps = {}
    for line in vcf_file:
        if line[0] != '#':
            split_line = line.split("\t")
            if len(split_line[3]) == 1 and len(split_line[4]) == 1 and split_line[4] != ".":
                tr_id = split_line[11]
                if tr_id not in _dict_snps:
                    _dict_snps[tr_id] = []
                _dict_snps[tr_id].append((split_line[2], split_line[0], int(split_line[1]), split_line[3], split_line[4], split_line[7]))
    vcf_file.close()
    return _dict_snps


def load_snp_table(data_path, version, GRCh, split=-1):
    snp_table = {}
    mk_data = open('{0}/{1}_{2}_estimates_snp{3}.out'.format(data_path, version, GRCh, "" if split == -1 else "_{0}".format(split)), 'r')
    mk_header = mk_data.readline().replace('\n', '').split('\t')
    for line in mk_data:
        line_split = line.replace('\n', '').split('\t')
        cds_id = line_split[0].split("_")[0]
        snp_table[cds_id] = dict(zip(mk_header[2:-2], [int(i) for i in line_split[2:-2]]))
        snp_table[cds_id][mk_header[1]] = line_split[1]
        assert len(line_split[1:-2]) == len(mk_header[1:-2])
        for index in [-1, -2]:
            if line_split[index] != "":
                snp_table[cds_id][mk_header[index]] = [int(i) for i in line_split[index].split(";")]
            else:
                snp_table[cds_id][mk_header[index]] = []
    mk_data.close()
    return snp_table

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

params_pb = [
    ("siteomega-predmutsel", "$\\omega_A = \\left< \\omega - \\omega_0 \\right>$"),
    ("predmutselfreeomega*(mutselfreeomega-1)", "$\\omega_A^* = \\left< \\omega_0^*  \\right>(\\omega^* - 1)$"),
    ("(siteomega-predmutsel)/siteomega", "$\\alpha = \\left< \\omega - \\omega_0 \\right> / \\left< \\omega \\right>$"),
    ("(mutselfreeomega-1)/mutselfreeomega", "$\\alpha^* = (\\omega^* - 1) / \\omega^* $"),
    ("predmutsel", "$\\left< \\omega_0 \\right>$"),
    ("predmutselfreeomega", "$\\left< \\omega_0^* \\right>$"),
    ("siteomega", "$\\left< \\omega \\right>$"),
    ("pos", "$\\left< P(\\omega > 1 ) \\right>$")]

columns = sorted(["siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega", "pos"])


def str_to_table(table, label):
    return eval(label, {f: table[f] for f in columns})


def load_pb_table(data_path):
    dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
    return np.loadtxt("{0}/79_GRCh38_estimates_pb.out".format(data_path), dtype=dtype, skiprows=1)

