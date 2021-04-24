#!/usr/bin/env python3
from collections import defaultdict
import gzip
import os
from lxml import etree
import pandas as pd
import numpy as np
from math import floor
from Bio.Phylo.PAML import yn00
from Bio import SeqIO, Seq

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


def build_divergence_dico(folder, gene_level=True, ci="0.025"):
    print('Loading divergence results.')
    ci = "0.0025" if gene_level else "0.05"
    dico_omega_0, dico_omega = {}, {}
    for file in sorted(os.listdir(folder)):
        sitemutsel_path = "{0}/{1}/sitemutsel_1.run.ci{2}.tsv".format(folder, file, ci)
        siteomega_path = "{0}/{1}/siteomega_1.run.ci{2}.tsv".format(folder, file, ci)
        if not os.path.isfile(siteomega_path) or not os.path.isfile(sitemutsel_path):
            continue

        engs = file[:-3]
        if gene_level:
            dico_omega_0[engs] = pd.read_csv(sitemutsel_path, sep="\t", nrows=1).values[:, 1:][0]
            dico_omega[engs] = pd.read_csv(siteomega_path, sep="\t", nrows=1).values[:, 1:][0]
        else:
            dico_omega_0[engs] = pd.read_csv(sitemutsel_path, sep="\t", skiprows=1).values[:, 1:]
            dico_omega[engs] = pd.read_csv(siteomega_path, sep="\t", skiprows=1).values[:, 1:]
            '''
            assert ((dico_omega_0[engs][:, 0] <= dico_omega_0[engs][:, 1]).all())
            assert ((dico_omega_0[engs][:, 1] <= dico_omega_0[engs][:, 2]).all())
            assert ((dico_omega[engs][:, 1] <= dico_omega[engs][:, 2]).all())
            assert ((dico_omega[engs][:, 1] <= dico_omega[engs][:, 2]).all())
            '''
    print('Divergence results loaded.')
    return dico_omega_0, dico_omega


def split_outliers(dico_omega_0, dico_omega, gene_level=True, filter_set=False):
    strongly_adaptive_dico, adaptive_dico, epistasis_dico, nearly_neutral_dico, unclassify_dico = {}, {}, {}, {}, {}

    for ensg, omega in dico_omega.items():
        if filter_set and ensg not in filter_set: continue
        omega_0 = dico_omega_0[ensg]
        if gene_level:
            if (omega[2] - omega[0] > 0.3) or (omega_0[2] - omega_0[0] > 0.3):
                continue
            if 1.0 < omega[0]:
                strongly_adaptive_dico[ensg] = None
            elif omega_0[2] < omega[0]:
                adaptive_dico[ensg] = None
            elif 0.05 < omega[0] and omega[2] < omega_0[0]:
                epistasis_dico[ensg] = None
            elif 0.05 <= omega[0] <= 1.0 and (
                    (omega_0[1] <= omega[1] <= omega_0[2] and omega[0] <= omega_0[1] <= omega[1]) or
                    (omega[1] <= omega_0[1] <= omega[2] and omega_0[0] <= omega[1] <= omega_0[1])):
                nearly_neutral_dico[ensg] = None
            else:
                unclassify_dico[ensg] = None
        else:
            strongly_adaptive_dico[ensg] = [i for i, w in enumerate(omega) if 1.0 < w[0] and w[1] < 1.4]
            adaptive_dico[ensg] = [i for i, w in enumerate(omega) if omega_0[i][2] < w[0] <= 1.0]
            epistasis_dico[ensg] = [i for i, w in enumerate(omega) if 0.05 < w[0] and w[2] < omega_0[i][0]]
            nearly_neutral_dico[ensg] = [i for i, w in enumerate(omega) if 0.05 <= w[0] <= 1.0 and (
                    (omega_0[i][1] <= w[1] <= omega_0[i][2] and w[0] <= omega_0[i][1] <= w[1]) or
                    (w[1] <= omega_0[i][1] <= w[2] and omega_0[i][0] <= w[1] <= omega_0[i][1]))]

            assert (len(set(nearly_neutral_dico[ensg]) & set(epistasis_dico[ensg])) == 0)
            assert (len(set(nearly_neutral_dico[ensg]) & set(adaptive_dico[ensg])) == 0)
            assert (len(set(nearly_neutral_dico[ensg]) & set(strongly_adaptive_dico[ensg])) == 0)
            assert (len(set(adaptive_dico[ensg]) & set(strongly_adaptive_dico[ensg])) == 0)
            assert (len(set(adaptive_dico[ensg]) & set(epistasis_dico[ensg])) == 0)

            classified = set(strongly_adaptive_dico[ensg]).union(set(adaptive_dico[ensg])).union(
                set(epistasis_dico[ensg])).union(set(nearly_neutral_dico[ensg]))
            unclassify_dico[ensg] = [i for i, w in enumerate(omega) if 0.05 <= w[0] and i not in classified]
    return strongly_adaptive_dico, adaptive_dico, epistasis_dico, nearly_neutral_dico, unclassify_dico


def filtered_table_omega(dico_omega, dico_subset, gene_level=True):
    output = []
    for ensg in dico_subset:
        ens_omega = dico_omega[ensg]
        if gene_level:
            output.append(ens_omega)
        else:
            output.extend([ens_omega[pos] for pos in dico_subset[ensg]])
    return np.array(output, dtype=np.float)


def table_omega(dico_omega, gene_level=True):
    output = []
    for values in dico_omega.values():
        if gene_level:
            output.append(values)
        else:
            output.extend(values)
    return np.array(output, dtype=np.float)


def snp_data_frame(vcf_path):
    print("Loading file " + vcf_path)
    fixed_poly = dict()
    sample_size_set = set()
    snp_table, header = {"ENSG": [], "POS": [], "TYPE": [], "COUNT": []}, {}
    vcf_file = gzip.open(vcf_path, 'rt')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                header = {k: i for i, k in enumerate(vcf_line.strip().split("\t"))}
            continue

        line_list = vcf_line.strip().split("\t")
        ensg = line_list[header["ENSG"]]
        if ensg == "None": continue

        if line_list[header["CHR"]] in ["X", "Y", "MT"]: continue

        genotypes = [s.split(":")[0].count("1") for s in line_list[header["FORMAT"] + 1:header["CHR"]] if
                     ("|" in s) or ("/" in s)]
        count = sum(genotypes)
        n = len(genotypes) * 2
        if count == 0: continue

        nuc_pos = int(line_list[header["ENSG_POS"]])
        codon_pos = int(nuc_pos / 3)
        if count == n:
            if ensg not in fixed_poly: fixed_poly[ensg] = {}
            ref, alt = line_list[header["REF"]], line_list[header["ALT"]]
            fixed_poly[ensg][codon_pos] = nuc_pos % 3, ref, alt
            continue

        snp_table["COUNT"].append(count)
        snp_table["ENSG"].append(ensg)
        snp_table["POS"].append(codon_pos)
        snp_table["TYPE"].append(line_list[header["SNP_TYPE"]])
        sample_size_set.add(n)

    print(vcf_path + " loaded.")
    assert len(sample_size_set) == 1
    n = sample_size_set.pop()
    return pd.DataFrame(snp_table).groupby("ENSG"), fixed_poly, n


def load_alignments(ensg_list, sp_1, sp_2, ali_folder):
    ali_dico = dict()

    for ensg in ensg_list:
        sister_focal_seqs = dict()
        for f in SeqIO.parse(open(ali_folder + "/" + ensg + "_NT.fasta", 'r'), 'fasta'):
            if f.id not in [sp_1, sp_2]: continue
            sister_focal_seqs[f.id] = f.seq

        if sp_1 in sister_focal_seqs and sp_2 in sister_focal_seqs:
            ali_dico[ensg] = sister_focal_seqs

    return ali_dico


def filter_positions(sister_focal_seqs, sp_focal, fixed_poly, gene_level, positions):
    seq_dico = dict()

    for sp_id, seq in sister_focal_seqs.items():

        s_list = list()
        for pos in (range(len(seq) // 3) if gene_level else positions):
            codon = str(seq[pos * 3:pos * 3 + 3])

            if pos in fixed_poly and sp_id == sp_focal:
                frame, ref, alt = fixed_poly[pos]
                if codon[frame] != ref:
                    assert codon[frame] == complement[ref]
                    alt = complement[alt]

                codon = codon[:frame] + alt + codon[frame + 1:]

            if codontable[codon] == "X": codon = "---"
            s_list.append(codon)

        seq_dico[sp_id] = "".join(s_list)

    return seq_dico


def run_yn00(seq1, seq2, tmp_path, filepath):
    SeqIO.write([seq1, seq2], filepath + ".fasta", "fasta")
    yn = yn00.Yn00(alignment=filepath + ".fasta", working_dir=tmp_path, out_file=filepath + "_yn00.out")

    try:
        res = yn.run()[seq1.id][seq2.id]['YN00']
        dn = int(round(res['N'] * res['dN']))
        ds = int(round(res['S'] * res['dS']))
        return res['N'], dn, res['S'], ds
    except:
        return 0, 0, 0, 0


def dfe_alpha(filepath, df, n, ensg_dico_pos, gene_level, sp_1, sp_2, ali_dico, fixed_poly, tmp_path, dfe_path):
    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    s1, s2 = [], []

    for ensg in ensg_dico_pos:
        seq_dico = filter_positions(ali_dico[ensg], sp_1, fixed_poly[ensg] if ensg in fixed_poly else {},
                                    gene_level, ensg_dico_pos[ensg])

        s1.append(seq_dico[sp_1])
        s2.append(seq_dico[sp_2])

        if ensg not in df.groups: continue

        dff = df.get_group(ensg)
        if not gene_level:
            filt = dff["POS"].isin(ensg_dico_pos[ensg])
            dff = dff[filt]

        for row in dff.itertuples(index=False):
            if row.TYPE == "Syn":
                ps += 1
                sfs_s[row.COUNT] += 1
            elif row.TYPE == "NonSyn":
                pn += 1
                sfs_n[row.COUNT] += 1
            else:
                assert row.TYPE == "Stop"

    seq1 = SeqIO.SeqRecord(id=sp_1, name="", description="", seq=Seq.Seq("".join(s1)))
    seq2 = SeqIO.SeqRecord(id=sp_2, name="", description="", seq=Seq.Seq("".join(s2)))
    sites_n, dn, sites_s, ds = run_yn00(seq1, seq2, tmp_path, filepath.replace(".txt", ""))

    sfs_list = ["Summed", n]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        range_sfs = range(1, int(floor(n // 2)) + 1)
        assert len(range_sfs) * 2 == n
        sfs_list += [nbr_site] + [(sfs[i] + sfs[n - i]) if n - i != i else sfs[i] for i in range_sfs]

    sfs_list += [sites_n, dn, sites_s, ds]
    content = "{0}+{1} ({2} sites)".format(sp_1, sp_2, len(seq1.seq)) + "\n"
    content += "\t".join(map(str, sfs_list)) + "\n"

    sfs_file = open(filepath, 'w')
    sfs_file.write(content)
    sfs_file.close()

    out_filename = filepath.replace(".txt", ".csv")
    os.system("{0}  -in {1} -out {2} -model GammaExpo".format(dfe_path, filepath, out_filename))


def subsample_sites(cds_dico, nbr_sites, weights):
    dico_list = list(cds_dico)
    rand_dico = dict()
    interval = [0]
    for k, v in cds_dico.items():
        interval.append(interval[-1] + len(v))

    if nbr_sites > interval[-1]:
        nbr_sites = interval[-1]
    site_choices = sorted(np.random.choice(interval[-1], nbr_sites, replace=False, p=weights))

    insert = np.searchsorted(interval, site_choices, side='right') - 1

    for ins, site in zip(insert, site_choices):
        ensg = dico_list[ins]
        pos = cds_dico[ensg][site - interval[ins]]
        if ensg not in rand_dico: rand_dico[ensg] = []
        rand_dico[ensg].append(pos)
    return rand_dico


def subsample_genes(cds_dico, nbr_sites, weights):
    return {k: None for k in np.random.choice(list(cds_dico), nbr_sites, replace=False, p=weights)}


def bin_dataset(dico_omega_0, dico_omega, bins=10, gene_level=True, filter_set=False):
    bin_dicos = [dict() for _ in range(bins)]
    interval = np.linspace(0, 0.5, bins)
    for ensg, omega in dico_omega.items():
        if filter_set and ensg not in filter_set: continue
        if gene_level:
            omega_A = omega[1] - dico_omega_0[ensg][1]
            insert = np.searchsorted(interval, omega_A, side='right') - 1
            bin_dicos[insert][ensg] = None
        else:
            for pos, site_omega in enumerate(omega):
                omega_A = site_omega[1] - dico_omega_0[ensg][pos][1]
                insert = np.searchsorted(interval, omega_A, side='right') - 1
                if ensg not in bin_dicos[insert]:
                    bin_dicos[insert][ensg] = list()
                bin_dicos[insert][ensg].append(pos)
    return bin_dicos
