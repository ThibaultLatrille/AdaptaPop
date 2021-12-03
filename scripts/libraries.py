#!/usr/bin/env python3
from collections import defaultdict
import gzip
import os
import pandas as pd
import numpy as np
from math import floor
from Bio.Phylo.PAML import yn00
from Bio import SeqIO, Seq

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
nucleotides = list(sorted(complement.keys()))
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

grapes_cmd = "{0} -in {1}.dofe -out {2}.csv -model GammaExpo -no_div_data 1> {2}.out 2> {2}.err"
polyDFE_cmd = "{0} -d {1}.sfs -w -m C -e 1> {2}.out 2> {2}.err"


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
            return "NotInCds", "", ""
        frame = nt_position % 3
        codon_ref = fasta_seq[nt_position - frame:nt_position + 3 - frame]
        aa_ref = codontable[codon_ref]
        if self.strand == "-":
            ref_nuc = complement[ref_nuc]
            alt_nuc = complement[alt_nuc]
        codon_alt = codon_ref[:frame] + alt_nuc + codon_ref[frame + 1:]
        aa_alt = codontable[codon_alt]
        if aa_ref == '' or len(codon_ref) != 3 or aa_ref == '-' or aa_alt == '-':
            return "NotIdentified", codon_ref, codon_alt
        elif codon_ref[frame] != ref_nuc:
            return "RefDiff", codon_ref, codon_alt
        elif aa_ref == 'X':
            return "RefStop", codon_ref, codon_alt
        elif aa_alt == 'X':
            return "Stop", codon_ref, codon_alt
        elif aa_alt == aa_ref:
            return "Syn", codon_ref, codon_alt
        else:
            return "NonSyn", codon_ref, codon_alt

    def empty_exon(self):
        return sum([1 for exon_len in self.exons_length() if exon_len == 1]) > 0

    def exons_length(self):
        return [j - i + 1 for i, j in self.exons]

    def seq_length(self):
        return sum(self.exons_length())

    def not_coding(self):
        return self.seq_length() % 3 != 0

    def befile_lines(self, chr_prefix):
        lines = []
        for start, end in self.exons:
            lines.append("{0}{1}\t{2}\t{3}\t{4}\t0\t{5}\n".format(chr_prefix, self.chromosome, start, end, self.name,
                                                                  self.strand))
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
    from lxml import etree
    dico_trid = {}
    for file in os.listdir(xml_folder):
        root = etree.parse(xml_folder + "/" + file).getroot()
        for info in root.findall(".//infoCDS[@specy='{0}']".format(specie)):
            trid = str(info.find('ensidTr').text)
            assert trid not in dico_trid
            dico_trid[trid] = file.replace(".xml", "")
    print('TR_ID to ENSG conversion done.')
    return dico_trid


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


def build_divergence_dico(folder, ensg_list, gene_level=True):
    print('Loading divergence results.')
    ci = "0.0025" if gene_level else "0.05"
    dico_omega_0, dico_omega = {}, {}
    for engs in sorted(ensg_list):
        sitemutsel_path = "{0}/{1}_NT/sitemutsel_1.run.ci{2}.tsv".format(folder, engs, ci)
        siteomega_path = "{0}/{1}_NT/siteomega_1.run.ci{2}.tsv".format(folder, engs, ci)
        if not os.path.isfile(siteomega_path) or not os.path.isfile(sitemutsel_path):
            sitemutsel_path = sitemutsel_path.replace("_null_", "__")
            siteomega_path = siteomega_path.replace("_null_", "__")

        assert os.path.isfile(siteomega_path) and os.path.isfile(sitemutsel_path)

        if gene_level:
            dico_omega_0[engs] = pd.read_csv(sitemutsel_path, sep="\t", nrows=1).values[:, 1:][0]
            dico_omega[engs] = pd.read_csv(siteomega_path, sep="\t", nrows=1).values[:, 1:][0]
        else:
            dico_omega_0[engs] = pd.read_csv(sitemutsel_path, sep="\t", skiprows=1).values[:, 1:]
            dico_omega[engs] = pd.read_csv(siteomega_path, sep="\t", skiprows=1).values[:, 1:]
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
            if 1.0 < omega[1]:
                strongly_adaptive_dico[ensg] = None
            elif omega_0[2] < omega[0] and omega[1] <= 1.0:
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
            strongly_adaptive_dico[ensg] = [i for i, w in enumerate(omega) if 1.0 < w[1] < 5.0]
            adaptive_dico[ensg] = [i for i, w in enumerate(omega) if omega_0[i][2] < w[0] and w[1] <= 1.0]
            epistasis_dico[ensg] = [i for i, w in enumerate(omega) if 0.05 < w[0] and w[2] < omega_0[i][0]]
            nearly_neutral_dico[ensg] = [i for i, w in enumerate(omega) if 0.05 <= w[0] and w[1] <= 1.0 and (
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


def snp_data_frame(vcf_path, is_unfolded, subsample):
    fixed_poly = dict()
    sample_size_set = set()
    snp_table = []
    dtype = {"CHR": 'string', "REF": 'string', "ALT": 'string', "ANC": 'string', "COUNT": int,
             "SAMPLE_SIZE": int, "ENSG": 'string', "POS": int, "TYPE": 'string'}
    df_snps = pd.read_csv(vcf_path, compression="gzip", dtype=dtype)
    tot = {"FIXED": 0, "REF": 0, "ALT": 0, "UNDEFINED": 0, "ABS": 0}
    for index, row in df_snps.iterrows():
        if subsample > 0:
            if subsample > row["SAMPLE_SIZE"]:
                subsample = row["SAMPLE_SIZE"]
            row["COUNT"] = np.random.hypergeometric(row["COUNT"], row["SAMPLE_SIZE"] - row["COUNT"], subsample)
            row["SAMPLE_SIZE"] = subsample

        codon_pos = int(row["POS"] / 3)
        if is_unfolded:
            if row["ANC"] == row["REF"]:
                tot["REF"] += 1
            elif row["ANC"] == row["ALT"]:
                row["COUNT"] = row["SAMPLE_SIZE"] - row["COUNT"]
                if row["ENSG"] not in fixed_poly: fixed_poly[row["ENSG"]] = {}
                fixed_poly[row.ENSG][codon_pos] = row["POS"] % 3, row["REF"], row["ALT"]
                tot["ALT"] += 1
            else:
                tot["UNDEFINED"] += 1
                continue
        else:
            tot["REF"] += 1

        if row["COUNT"] == 0:
            tot["ABS"] += 1
            continue
        elif row["COUNT"] == row["SAMPLE_SIZE"]:
            if row["ENSG"] not in fixed_poly: fixed_poly[row["ENSG"]] = {}
            fixed_poly[row["ENSG"]][codon_pos] = row["POS"] % 3, row["REF"], row["ALT"]
            tot["FIXED"] += 1
            continue

        snp_table.append(row)
        sample_size_set.add(row["SAMPLE_SIZE"])
    tot_sum = sum(tot.values())
    print("{0} SNPs".format(tot_sum))
    for k, v in tot.items():
        print("{0}: {1} ({2:.2f}%)".format(k, v, 100 * v / tot_sum))
    assert len(sample_size_set) == 1
    return pd.DataFrame(snp_table).groupby("ENSG"), fixed_poly, sample_size_set.pop()


def load_alignments(ensg_list, sp_1, sp_2, ali_folder):
    ali_dico = dict()

    for ensg in ensg_list:
        sister_focal_seqs = dict()
        file_path = ali_folder + "/" + ensg + "_NT.fasta"
        if not os.path.exists(file_path): file_path = file_path.replace("_null_", "__")
        for f in SeqIO.parse(open(file_path, 'r'), 'fasta'):
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
                codon = codon[:frame] + alt + codon[frame + 1:]

            if codontable[codon] == "X" or codontable[codon] == "-":
                codon = "---"
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


def dfe_alpha(filepath, df, n, ensg_dico_pos, gene_level, sp_1, sp_2, ali_dico, fixed_poly, tmp_path, dfe_models,
              is_unfolded, error_f):
    if dfe_models is None:
        dfe_models = []
    out_list = [os.path.isfile(filepath + "_" + dfe_path.split("/")[-2] + (".csv" if "grapes" in dfe_path else ".out"))
                for dfe_path in dfe_models]
    if np.all(out_list) and len(dfe_models) > 0:
        print("Output files already exists, not computing and exiting")
        return True
    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    s1, s2 = [], []

    for ensg in ensg_dico_pos:
        seq_dico = filter_positions(ali_dico[ensg], sp_1, fixed_poly[ensg] if ensg in fixed_poly else {},
                                    gene_level, ensg_dico_pos[ensg])

        s1.append(seq_dico[sp_1])
        s2.append(seq_dico[sp_2])

        if ensg not in df.groups:
            continue

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

    if np.sum(sfs_n) < 10 or np.sum(sfs_s) < 10:
        error_f.write("ER0: not enough polymorphism. SFS are:\n" + " ".join(sfs_s) + "\n" + " ".join(sfs_n) + "\n")
        return False

    seq1 = SeqIO.SeqRecord(id=sp_1, name="", description="", seq=Seq.Seq("".join(s1)))
    seq2 = SeqIO.SeqRecord(id=sp_2, name="", description="", seq=Seq.Seq("".join(s2)))
    if seq1.seq == seq2.seq:
        error_f.write("ER1: identical sequences. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    sites_n, dn, sites_s, ds = run_yn00(seq1, seq2, tmp_path, filepath)

    if sites_n == dn == sites_s == ds == 0:
        print('error in Yn00')
        error_f.write("ER2: Yn00 failed. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    sfs_list = ["Summed", n]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        if is_unfolded:
            sfs_list += [nbr_site] + [sfs[i] for i in range(1, n)]
        else:
            range_sfs = range(1, int(floor(n // 2)) + 1)
            assert len(range_sfs) * 2 == n
            sfs_list += [nbr_site] + [(sfs[i] + sfs[n - i]) if n - i != i else sfs[i] for i in range_sfs]
    sfs_list += [sites_n, dn, sites_s, ds]

    dofe_file = open(filepath + ".dofe", 'w')
    dofe_file.write("{0}+{1} ({2} sites)".format(sp_1, sp_2, len(seq1.seq)) + "\n")
    if is_unfolded:
        dofe_file.write("#unfolded\n")
    dofe_file.write("\t".join(map(str, sfs_list)) + "\n")
    dofe_file.close()

    sfs_s_str = " ".join([str(sfs_s[i]) for i in range(1, n)]) + "\t{0}\t{1}\t{0}\n".format(sites_s, ds)
    sfs_n_str = " ".join([str(sfs_n[i]) for i in range(1, n)]) + "\t{0}\t{1}\t{0}\n".format(sites_n, dn)
    sfs_file = open(filepath + ".sfs", 'w')
    sfs_file.write("#{0}+{1}".format(sp_1, sp_2) + "\n")
    sfs_file.write("1 1 {0}".format(n) + "\n")
    sfs_file.write(sfs_s_str)
    sfs_file.write(sfs_n_str)
    sfs_file.close()

    for dfe_path in dfe_models:
        out = filepath + "_" + dfe_path.split("/")[-2]
        if "grapes" in dfe_path:
            os.system(grapes_cmd.format(dfe_path, filepath, out))
        elif "polyDFE" in dfe_path and is_unfolded:
            os.system(polyDFE_cmd.format(dfe_path, filepath, out))
    os.system("gzip --force {0}.fasta".format(filepath))
    return True


def subsample_sites(cds_dico, nbr_sites, weights, replace=False):
    dico_list = list(cds_dico)
    rand_dico = dict()
    interval = [0]
    for k, v in cds_dico.items():
        interval.append(interval[-1] + len(v))

    if nbr_sites > interval[-1]:
        nbr_sites = interval[-1]
    site_choices = sorted(np.random.choice(interval[-1], nbr_sites, replace=replace, p=weights))

    insert = np.searchsorted(interval, site_choices, side='right') - 1

    for ins, site in zip(insert, site_choices):
        ensg = dico_list[ins]
        pos = cds_dico[ensg][site - interval[ins]]
        if ensg not in rand_dico: rand_dico[ensg] = []
        rand_dico[ensg].append(pos)
    return rand_dico


def subsample_genes(cds_dico, nbr_genes, weights, replace=False):
    return {k: None for k in np.random.choice(list(cds_dico), nbr_genes, replace=replace, p=weights)}


GREEN = "#8FB03E"
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"


def format_pop(t):
    if " " in t:
        return "".join([s[0] for s in t.split(" ")])
    else:
        return t


def sp_to_color(specie):
    if "Homo" in specie:
        return BLUE
    elif "Chloro" in specie:
        return "black"
    elif "Ovis" in specie:
        return YELLOW
    elif "Bos" in specie:
        return LIGHTGREEN
    else:
        return "grey"


def sp_sorted(pop, sp):
    out = sp + "_" + pop
    if "Homo" in out:
        return "Z" + out
    elif "Chloro" in out:
        return "Y" + out
    elif "Ovis" in out:
        return "X" + out
    elif "Bos" in out:
        return "W" + out
    else:
        return out
