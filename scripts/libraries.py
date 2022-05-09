#!/usr/bin/env python3
from collections import defaultdict, namedtuple
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

confidence_interval = namedtuple('confidence_interval', ['low', 'mean', 'up'])
sfs_weight = {"watterson": lambda i, n: 1.0 / i, "tajima": lambda i, n: n - i, "fay_wu": lambda i, n: i}


def theta(sfs_epsilon, daf_n, weight_method):
    sfs_theta = sfs_epsilon * np.array(range(1, daf_n))
    weights = np.array([sfs_weight[weight_method](i, daf_n) for i in range(1, daf_n)])
    return sum(sfs_theta * weights) / sum(weights)


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

    def overlapping(self):
        for i in range(len(self.exons) - 1):
            if not self.exons[i][1] < self.exons[i + 1][0]:
                return True
        return False

    def befile_lines(self, chr_prefix):
        lines = []
        for start, end in self.exons:
            lines.append("{0}{1}\t{2}\t{3}\t{4}\t0\t{5}\n".format(chr_prefix, self.chromosome, start, end, self.name,
                                                                  self.strand))
        return lines


def build_dict_cds(data_path, file_name, chr2acc=""):
    dico_chr2acc = dict()
    if chr2acc != "":
        chr2acc_file = open("{0}/{1}".format(data_path, chr2acc), 'rt')
        chr2acc_file.readline()
        for line in chr2acc_file:
            splitline = line.strip().split("\t")
            dico_chr2acc[splitline[1]] = splitline[0]

    print('Loading GTF file...')
    gtf_file = gzip.open("{0}/{1}".format(data_path, file_name), 'rt')
    dico_cds, pr2tr_id = dict(), dict()
    nf_tr_id = set()
    for gtf_line in gtf_file:
        if gtf_line.startswith('#'):
            continue

        chromosome, source, feature, start, end, score, strand, frame, comments = gtf_line.replace('\n', '').split('\t')
        if feature != 'CDS':
            continue

        if chromosome in dico_chr2acc:
            chromosome = dico_chr2acc[chromosome]
        transcript_find = comments.find('transcript_id')
        # comments.find('CCDS') != -1 and comments.find('ccds_id') != -1
        if transcript_find != -1:
            tr_id = comments[transcript_find + 15:].split("\"")[0]
            if comments.find('cds_start_NF') != -1 or comments.find('cds_end_NF') != -1:
                if tr_id not in nf_tr_id:
                    nf_tr_id.add(tr_id)
            if tr_id not in dico_cds:
                dico_cds[tr_id] = Cds(chromosome, strand, tr_id)
            dico_cds[tr_id].add_exon(start, end)
            protein_find = comments.find('protein_id')
            if protein_find != -1:
                pr_id = comments[protein_find + 12:].split("\"")[0]
                if pr_id in pr2tr_id:
                    assert pr2tr_id[pr_id] == tr_id
                pr2tr_id[pr_id] = tr_id
    gtf_file.close()
    print('GTF file loaded.')
    return dico_cds, nf_tr_id, pr2tr_id


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


def build_divergence_dico(folder, ensg_list, gene_level=True, pp="0.025"):
    print('Loading divergence results.')
    assert pp in ["0.025", "0.0025"]
    dico_omega_0, dico_omega = {}, {}
    for engs in sorted(ensg_list):
        sitemutsel_path = "{0}/{1}_NT/sitemutsel_1.run.ci{2}.tsv".format(folder, engs, pp)
        siteomega_path = "{0}/{1}_NT/siteomega_1.run.ci{2}.tsv".format(folder, engs, pp)
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


def is_adaptive(omega, omega_0, method):
    return (method == "Classical" and (1.0 < omega.low)) or (
            method == "MutSel" and (omega_0.up < omega.low)) or (
                   method == "MutSelExclu" and (omega_0.up < omega.low and omega.mean <= 1.0))


def split_outliers(dico_omega_0, dico_omega, gene_level=True, filter_set=False, method="MutSel"):
    assert method in ["MutSel", "Classical", "MutSelExclu"]
    adaptive_dico, nearly_neutral_dico, unclassify_dico = {}, {}, {}

    for ensg in dico_omega:
        if filter_set and (ensg not in filter_set):
            continue

        if gene_level:
            omega = confidence_interval(*dico_omega[ensg])
            omega_0 = confidence_interval(*dico_omega_0[ensg])
            if (omega.up - omega.low > 0.3) or (omega_0.up - omega_0.low > 0.3):
                unclassify_dico[ensg] = None
            elif is_adaptive(omega, omega_0, method):
                adaptive_dico[ensg] = None
            elif (0.05 <= omega.low and omega.mean <= 1.0) and (
                    (omega_0.mean <= omega.mean <= omega_0.up and omega.low <= omega_0.mean <= omega.mean) or
                    (omega.mean <= omega_0.mean <= omega.up and omega_0.low <= omega.mean <= omega_0.mean)):
                nearly_neutral_dico[ensg] = None
            else:
                unclassify_dico[ensg] = None
        else:
            w_list = [confidence_interval(*i) for i in dico_omega[ensg]]
            w_0_list = [confidence_interval(*i) for i in dico_omega_0[ensg]]
            adaptive_dico[ensg] = [i for i, w in enumerate(w_list) if is_adaptive(w, w_0_list[i], method)]
            nearly_neutral_dico[ensg] = [i for i, w in enumerate(w_list) if 0.05 <= w.low and w.mean <= 1.0 and (
                    (w_0_list[i].mean <= w.mean <= w_0_list[i].up and w.low <= w_0_list[i].mean <= w.mean) or
                    (w.mean <= w_0_list[i].mean <= w.up and w_0_list[i].low <= w.mean <= w_0_list[i].mean))]

            assert (len(set(nearly_neutral_dico[ensg]) & set(adaptive_dico[ensg])) == 0)
            classified = set(adaptive_dico[ensg]).union(set(nearly_neutral_dico[ensg]))
            unclassify_dico[ensg] = [i for i in range(len(w_list)) if i not in classified]

    return adaptive_dico, nearly_neutral_dico, unclassify_dico


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


def snp_data_frame(vcf_path, polarize_snps, remove_fixed=True):
    dtype = {"CHR": 'string', "REF": 'string', "ALT": 'string', "ANC": 'string', "ANC_PROBA": float, "COUNT": int,
             "SAMPLE_SIZE": int, "ENSG": 'string', "CODON_POS": int, "NUC_POS": int, "TYPE": 'string'}
    df_snps = pd.read_csv(vcf_path, compression="gzip", dtype=dtype)
    fixed_poly = defaultdict(dict)
    if polarize_snps:
        sample_size_set = set()
        snp_table = []
        tot = defaultdict(lambda: 0)
        for index, row in df_snps.iterrows():
            is_polarized = row["ANC"] == row["REF"]
            if is_polarized:
                tot["REF=ANC"] += 1
                anc, der = row["REF"], row["ALT"]
                count_polarized = row["COUNT"]
            elif row["ANC"] == row["ALT"]:
                # Flip the count and the alt is fixed
                count_polarized = row["SAMPLE_SIZE"] - row["COUNT"]
                anc, der = row["ALT"], row["REF"]
                fixed_poly[row["ENSG"]][row["CODON_POS"]] = row["NUC_POS"] % 3, row["REF"], anc
                tot["ALT=ANC"] += 1
            else:
                tot["ANC_UNDEFINED"] += 1
                continue

            if count_polarized == 0:
                tot["ABS"] += 1
                fixed_poly[row["ENSG"]][row["CODON_POS"]] = row["NUC_POS"] % 3, row["REF"], anc
                if remove_fixed:
                    continue
            elif count_polarized == row["SAMPLE_SIZE"]:
                fixed_poly[row["ENSG"]][row["CODON_POS"]] = row["NUC_POS"] % 3, row["REF"], der
                tot["FIXED"] += 1
                if remove_fixed:
                    continue

            row["DER"], row["POLARIZED"], row["COUNT_POLARIZED"] = der, is_polarized, count_polarized
            snp_table.append(row)
            sample_size_set.add(row["SAMPLE_SIZE"])

        tot_sum = sum(tot.values())
        print("{0} SNPs".format(tot_sum))
        for k, v in tot.items():
            print("{0}: {1} ({2:.2f}%)".format(k, v, 100 * v / tot_sum))
        df_snps = pd.DataFrame(snp_table)
    else:
        sample_size_set = set(df_snps["SAMPLE_SIZE"])

    assert len(sample_size_set) == 1
    return df_snps.groupby("ENSG"), fixed_poly, sample_size_set.pop()


def load_alignments(ensg_list, sp_1, sp_2, ali_folder):
    ali_dico = dict()

    for ensg in ensg_list:
        sister_focal_seqs = dict()
        file_path = ali_folder + "/" + ensg + "_NT.fasta"
        if not os.path.exists(file_path):
            file_path = file_path.replace("_null_", "__")
        for f in SeqIO.parse(open(file_path, 'r'), 'fasta'):
            if f.id not in [sp_1, sp_2]:
                continue
            sister_focal_seqs[f.id] = f.seq

        if sp_1 in sister_focal_seqs and sp_2 in sister_focal_seqs:
            ali_dico[ensg] = sister_focal_seqs

    return ali_dico


def filter_positions(sister_focal_seqs, sp_focal, fixed_poly, gene_level, positions):
    seq_dico = dict()

    for sp_id, seq in sister_focal_seqs.items():

        s_list = list()
        for codon_pos in (range(len(seq) // 3) if gene_level else positions):
            codon = str(seq[codon_pos * 3:codon_pos * 3 + 3])

            if codon_pos in fixed_poly and sp_id == sp_focal:
                frame, ref, alt = fixed_poly[codon_pos]
                assert codon[frame] == ref
                codon = codon[:frame] + alt + codon[frame + 1:]

            if codontable[codon] == "X" or codontable[codon] == "-":
                codon = "---"
            s_list.append(codon)

        seq_dico[sp_id] = "".join(s_list)

    return seq_dico


def run_yn00(seq1, seq2, tmp_path, filepath):
    SeqIO.write([seq1, seq2], f"{filepath}.fasta", "fasta")
    yn = yn00.Yn00(alignment=f"{filepath}.fasta", working_dir=tmp_path, out_file=f"{filepath}_yn00.out")

    try:
        res = yn.run()[seq1.id][seq2.id]['YN00']
        dn = int(round(res['N'] * res['dN']))
        ds = int(round(res['S'] * res['dS']))
        return res['N'], dn, res['S'], ds
    except:
        return 0, 0, 0, 0


def write_dofe(sfs_syn, sfs_non_syn, l_non_syn, d_non_syn, l_syn, d_syn, k, filepath, sp_focal, sp_sister, is_unfolded,
               L):
    sfs_list = [k]
    for sfs, nbr_site in [(sfs_non_syn, l_non_syn), (sfs_syn, l_syn)]:
        if is_unfolded:
            sfs_list += [nbr_site] + [sfs[i] for i in range(1, k)]
        else:
            range_sfs = range(1, int(floor(k // 2)) + 1)
            assert len(range_sfs) * 2 == k
            sfs_list += [nbr_site] + [(sfs[i] + sfs[k - i]) if k - i != i else sfs[i] for i in range_sfs]
    sfs_list += [l_non_syn, d_non_syn, l_syn, d_syn]

    dofe_file = open(filepath + ".dofe", 'w')
    dofe_file.write(f"{sp_focal}+{sp_sister} ({int(L)} sites)\n")
    if is_unfolded:
        dofe_file.write("#unfolded\n")
    dofe_file.write("Summed\t" + "\t".join(map(lambda i: str(int(i)), sfs_list)) + "\n")
    dofe_file.close()


def grapes_cmd(dfe_path, filepath, out):
    return f"{dfe_path} -in {filepath}.dofe -out {out}.csv 1> {out}.out 2> {out}.err"


def write_sfs(sfs_syn, sfs_non_syn, l_non_syn, d_non_syn, l_syn, d_syn, k, filepath, sp_focal, sp_sister, div=True):
    sfs_syn_str = " ".join([str(int(sfs_syn[i])) for i in range(1, k)]) + f"\t{int(l_syn)}"
    sfs_non_syn_str = " ".join([str(int(sfs_non_syn[i])) for i in range(1, k)]) + f"\t{int(l_non_syn)}"
    if div:
        sfs_syn_str += f"\t{int(d_syn)}\t{int(l_syn)}"
        sfs_non_syn_str += f"\t{int(d_non_syn)}\t{int(l_non_syn)}"
    sfs_file = open(filepath + ".sfs", 'w')
    sfs_file.write(f"#{sp_focal}+{sp_sister}\n")
    sfs_file.write("1 1 {0}".format(k) + "\n")
    sfs_file.write(sfs_syn_str + "\n")
    sfs_file.write(sfs_non_syn_str + "\n")
    sfs_file.close()


def polyDFE_cmd(dfe_path, filepath, out):
    return f"{dfe_path} -d {filepath}.sfs 1> {out}.out 2> {out}.err"


def dfe_alpha(filepath, df, k_tot, k, ensg_dico_pos, gene_level, sp_focal, sp_sister, ali_dico, fixed_poly, tmp_path,
              dfe_models, is_unfolded, error_f):
    if dfe_models is None:
        dfe_models = []

    p_non_syn, p_syn = 0, 0
    sfs_non_syn, sfs_syn = np.zeros(k, dtype=int), np.zeros(k, dtype=int)
    s1, s2 = [], []

    for ensg in ensg_dico_pos:
        assert ensg in df.groups
        ensg_fixedpoly = dict(fixed_poly[ensg])

        dff = df.get_group(ensg)
        if not gene_level:
            filt = dff["CODON_POS"].isin(ensg_dico_pos[ensg])
            dff = dff[filt]
            assert id(dff) != id(df.get_group(ensg))

        for row in dff.itertuples(index=False):
            if row.TYPE == "Stop":
                continue

            if k_tot > k:
                daf = np.random.hypergeometric(row.COUNT, k_tot - row.COUNT, k)
                if is_unfolded and ((daf == 0) or (daf == k)):
                    if row.REF == row.ANC:
                        anc, der = row.REF, row.ALT
                    else:
                        anc, der = row.ALT, row.REF
                    if daf == 0:
                        ensg_fixedpoly[row.CODON_POS] = row.NUC_POS % 3, row.REF, anc
                        continue
                    elif daf == k:
                        ensg_fixedpoly[row.CODON_POS] = row.NUC_POS % 3, row.REF, der
                        continue
                elif (daf == 0) or (daf == k):
                    continue
            else:
                daf = row.COUNT

            if row.TYPE == "Syn":
                p_syn += 1
                sfs_syn[daf] += 1
            else:
                assert row.TYPE == "NonSyn"
                p_non_syn += 1
                sfs_non_syn[daf] += 1

        seq_dico = filter_positions(ali_dico[ensg], sp_focal, ensg_fixedpoly, gene_level, ensg_dico_pos[ensg])
        s1.append(seq_dico[sp_focal])
        s2.append(seq_dico[sp_sister])

    assert sfs_non_syn[0] == 0
    assert sfs_syn[0] == 0

    if np.sum(sfs_non_syn) < 10 or np.sum(sfs_syn) < 10:
        error_f.write(
            "ER0: not enough polymorphism. SFS are:\n" + " ".join(sfs_syn) + "\n" + " ".join(sfs_non_syn) + "\n")
        return False

    seq1 = SeqIO.SeqRecord(id=sp_focal, name="", description="", seq=Seq.Seq("".join(s1)))
    seq2 = SeqIO.SeqRecord(id=sp_sister, name="", description="", seq=Seq.Seq("".join(s2)))
    if seq1.seq == seq2.seq:
        error_f.write("ER1: identical sequences. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    l_non_syn, d_non_syn, l_syn, d_syn = run_yn00(seq1, seq2, tmp_path, filepath)

    if l_non_syn == d_non_syn == l_syn == d_syn == 0:
        print('error in Yn00')
        error_f.write("ER2: Yn00 failed. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    write_dofe(sfs_syn, sfs_non_syn, l_non_syn, d_non_syn, l_syn, d_syn, k,
               filepath, sp_focal, sp_sister, is_unfolded, len(seq1.seq))
    write_sfs(sfs_syn, sfs_non_syn, l_non_syn, d_non_syn, l_syn, d_syn, k, filepath, sp_focal, sp_sister)

    for dfe_path in dfe_models:
        model = dfe_path[:dfe_path.find(':')]
        dfe_cmd = dfe_path[dfe_path.find(':') + 1:]
        out = f"{filepath}_{model}"
        if "grapes" == model:
            cmd = grapes_cmd(dfe_cmd, filepath, out)
        else:
            assert "polyDFE" == model and is_unfolded
            cmd = polyDFE_cmd(dfe_cmd, filepath, out)
        print(cmd)
        os.system(cmd)
    os.system("gzip --force {0}.fasta".format(filepath))
    return True


def subsample_sites(cds_dico, nbr_sites, weights, replace=False):
    dico_list = list(cds_dico)
    rand_dico = defaultdict(list)
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
        rand_dico[ensg].append(pos)
    return rand_dico


def subsample_genes(cds_dico, nbr_genes, weights, replace=False):
    return {k: None for k in np.random.choice(list(cds_dico), nbr_genes, replace=replace, p=weights)}
