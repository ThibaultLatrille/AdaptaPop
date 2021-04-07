#!/usr/bin/env python3
from collections import defaultdict
import gzip
import os
from lxml import etree
import pandas as pd
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import numpy as np

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
nucleotides = list(complement.keys())
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
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W', '---': '-'})


def translate(codon_seq, gap=False):
    prot = []

    for n in range(0, len(codon_seq), 3):
        codon = codon_seq[n:n + 3]
        if codon in codontable:
            residue = codontable[codon]
        else:
            residue = "X"

        if gap:
            residue = " " + residue + " "
        prot.append(residue)

    return "".join(prot)


class Cds(object):
    def __init__(self, chromosome, strand, name):
        """
        :chromosome : (String) Chromosome number (can also be X/Y or Z/W).
        :strand : (String) Strand on which the CDS is encoded.
        :name : (String) Name of the CDS.
        :exons : (List of 2-tuple) List of exons. Each exon is defined by a tuple (start, end),
                 where 'start' and 'end' are in absolute position in the chromosome.
        """
        self.chromosome = chromosome
        assert (strand == "+" or strand == "-")
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        if int(start_exon) <= int(end_exon):
            if self.strand == "+":
                self.exons.append((int(start_exon), int(end_exon)))
            else:
                self.exons.insert(0, (int(start_exon), int(end_exon)))
            for i in range(len(self.exons) - 1):
                if not self.exons[i][1] < self.exons[i + 1][0]:
                    print("At least one exon is overlapping with an other")

    def nt_position(self, position):
        """
        :param position: (Integer) Nucleotide position in the chromosome.
        :return: (Integer) Nucleotide position relative to the CDS.
        """
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
        return nt

    def snp_type(self, fasta_seq, position, ref_nuc, alt_nuc):
        """
        Return the SNP type given the position, reference and alternative nucleotide.
        This method also requires the fasta sequence as input.
        :param fasta_seq: (String) The fasta sequence.
        :param position: (Integer) Absolute position of the SNP in the chromosome.
        :param ref_nuc: (Character) The reference nucleotide.
        :param alt_nuc: (Character) The absolute nucleotide.
        :return: (String) The classification of the SNP, can be either "NotInCds",
                          "RefDiff", "NotIdentified", "RefStop", "Stop", "Syn" or "NonSyn".
        """
        nt_position = self.nt_position(position)
        if nt_position == -1:
            return "NotInCds"
        frame = nt_position % 3
        codon_ref = fasta_seq[nt_position - frame:nt_position + 3 - frame]
        aa_ref = codontable[codon_ref]
        if self.strand == "-":
            ref_nuc = complement[ref_nuc]
            alt_nuc = complement[alt_nuc]
        codon_alt = codon_ref[:frame] + alt_nuc + codon_ref[frame + 1:]
        aa_alt = codontable[codon_alt]
        if aa_ref == '' or len(codon_ref) != 3 or aa_ref == '-' or aa_alt == '-':
            return "NotIdentified"
        elif codon_ref[frame] != ref_nuc:
            return "RefDiff"
        elif aa_ref == 'X':
            return "RefStop"
        elif aa_alt == 'X':
            return "Stop"
        elif aa_alt == aa_ref:
            return "Syn"
        else:
            return "NonSyn"

    def empty_exon(self):
        return sum([1 for exon_len in self.exons_length() if exon_len == 1]) > 0

    def exons_length(self):
        return [j - i + 1 for i, j in self.exons]

    def seq_length(self):
        return sum(self.exons_length())

    def not_coding(self):
        return self.seq_length() % 3 != 0

    def befile_lines(self):
        lines = []
        for start, end in self.exons:
            lines.append("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(self.chromosome, start, end, self.name, self.strand))
        return lines


def build_dict_cds(data_path, file_name):
    print('Loading GTF file...')
    gtf_file = gzip.open("{0}/{1}".format(data_path, file_name), 'rt')
    dico_cds = dict()
    nf_tr_id = set()
    for gtf_line in gtf_file:
        if gtf_line.startswith('#'):
            continue

        seqname, source, feature, start, end, score, strand, frame, comments = gtf_line.replace('\n', '').split('\t')
        if feature != 'CDS':
            continue

        transcript_find = comments.find('transcript_id')
        # comments.find('CCDS') != -1 and comments.find('ccds_id') != -1
        if transcript_find != -1:
            tr_id = comments[transcript_find + 15:].split("\"")[0]
            if comments.find('cds_start_NF') != -1 or comments.find('cds_end_NF') != -1:
                if tr_id not in nf_tr_id:
                    nf_tr_id.add(tr_id)
            if tr_id not in dico_cds:
                dico_cds[tr_id] = Cds(seqname, strand, tr_id)
            dico_cds[tr_id].add_exon(start, end)
    gtf_file.close()
    print('GTF file loaded.')
    return dico_cds, nf_tr_id


def build_dict_trID(xml_folder, specie):
    print('Converting TR_ID to ENSG.')
    dico_trid = {}
    for file in os.listdir(xml_folder):
        root = etree.parse(xml_folder + "/" + file).getroot()
        for info in root.findall(".//infoCDS[@specy='{0}']".format(specie)):
            trid = str(info.find('ensidTr').text)
            assert trid not in dico_trid
            dico_trid[trid] = file.replace(".xml", "")
    print('TR_ID to ENSG conversion done.')
    return dico_trid


def ontology_table(xml_folder):
    print('Finding CDS ontologies.')
    go_id2name, go_id2cds_list = {}, {}
    all_go_set = set()
    for file in os.listdir(xml_folder):
        root = etree.parse(xml_folder + "/" + file).getroot()
        for annot in root.find('goAnnot').findall("annot"):
            go_id = annot.find('goId').text
            go_name = annot.find('goName').text
            if go_id not in go_id2name: go_id2name[go_id] = go_name.replace('"', '')
            if go_id not in go_id2cds_list: go_id2cds_list[go_id] = set()
            ensg = file.replace(".xml", "")
            go_id2cds_list[go_id].add(ensg)
            all_go_set.add(ensg)
    print('CDS ontologies found.')
    return go_id2cds_list, go_id2name, all_go_set


def tex_f(x):
    if 0.001 < abs(x) < 10:
        return "{:6.3f}".format(x)
    elif 10 <= abs(x) < 10000:
        return "{:6.1f}".format(x)
    else:
        s = "{:6.2g}".format(x)
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + 'e$^{' + str(int(exp)) + '}$'
            s = " " * (5 + 6 - len(s)) + s
        return s


def build_divergence_dico(folder, gene_level=True):
    print('Loading divergence results.')
    dico_omega_0, dico_omega = {}, {}
    for file in sorted(os.listdir(folder)):
        sitemutsel_path = folder + "/" + file + "/sitemutsel_1.run.omegappgt1.000000"
        siteomega_path = folder + "/" + file + "/siteomega_1.run.omegappgt1.000000"
        if not os.path.isfile(siteomega_path) or not os.path.isfile(sitemutsel_path):
            continue

        engs = file[:-3]
        site_omega_0 = pd.read_csv(sitemutsel_path, sep="\t")["omega_0"].values
        site_omega = pd.read_csv(siteomega_path, sep="\t")["omega"].values
        if gene_level:
            dico_omega_0[engs] = np.mean(site_omega_0)
            dico_omega[engs] = np.mean(site_omega)
        else:
            assert len(site_omega_0) == len(site_omega)
            dico_omega_0[engs] = site_omega_0
            dico_omega[engs] = site_omega

    print('Divergence results loaded.')
    return dico_omega_0, dico_omega


def split_outliers(dico_omega_0, dico_omega, gene_level=True, filter_set=False):
    adaptive_dico, epistasis_dico, nearly_neutral_dico = {}, {}, {}

    for ensg, omega in dico_omega.items():
        if filter_set and ensg not in filter_set: continue
        omega_0 = dico_omega_0[ensg]
        if gene_level:
            if omega > omega_0 + 0.1:
                adaptive_dico[ensg] = None
            elif 0.05 < omega < omega_0 - 0.1:
                epistasis_dico[ensg] = None
            elif max(0.05, omega_0 - 0.1) <= omega <= omega_0 + 0.1:
                nearly_neutral_dico[ensg] = None
        else:
            adaptive_dico[ensg] = [i for i, v in enumerate(omega) if omega_0[i] + 0.1 < v < 1.0]
            epistasis_dico[ensg] = [i for i, v in enumerate(omega) if 0.05 < v < omega_0[i] - 0.1]
            nearly_neutral_dico[ensg] = [i for i, v in enumerate(omega) if
                                         max(0.05, omega_0[i] - 0.1) <= v <= omega_0[i] + 0.1]

    return adaptive_dico, epistasis_dico, nearly_neutral_dico


def filtered_table_omega(dico_omega, dico_subset, gene_level=True):
    output = []
    for ensg in dico_subset:
        ens_omega = dico_omega[ensg]
        if gene_level:
            output.append(ens_omega)
        else:
            output.extend([ens_omega[pos] for pos in dico_subset[ensg]])
    return output


def table_omega(dico_omega, gene_level=True):
    output = []
    for values in dico_omega.values():
        if gene_level:
            output.append(values)
        else:
            output.extend(values)
    return output
