#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.pairwise2 import align
from libraries import *


def most_common(lst):
    return max(set(lst), key=lst.count)


class Alignment(object):
    def __init__(self, cds_seq, omm_seq, ensg, file_err):
        self.error = ""
        self.ensg = ensg
        self.file_err = file_err
        if len(omm_seq) == 0:
            self.error = "SeqOMMEmpty"
            return

        self.omm_seq = omm_seq
        omm_seq_aligned = str(omm_seq.replace("---", ''))
        if len(omm_seq_aligned) % 3 != 0:
            self.error = "SeqOMMnotMultiple3"
            return

        if len(cds_seq) % 3 != 0:
            self.error = "SeqCDSnotMultiple3"
            return

        alns = align.globalmd(translate(cds_seq), translate(omm_seq_aligned), 10, -1000, -100, -10, -10, -1)
        if len(alns) == 0:
            self.error = "NoAlignment"
            return

        score = 0
        for aln in alns:
            seq1 = Alignment.gapsFromPeptide(str(aln[0]), cds_seq)
            seq2 = Alignment.gapsFromPeptide(str(aln[1]), omm_seq_aligned)
            tmp_score = Alignment.score_ali(seq1, seq2)
            if tmp_score >= score:
                self.cds_seq_aligned = seq1
                self.omm_seq_aligned = seq2

    def __repr__(self):
        return translate(self.cds_seq_aligned, True) + '\n' + self.cds_seq_aligned + '\n' + \
               self.omm_seq_aligned + '\n' + translate(self.omm_seq_aligned, True) + '\n' + self.omm_seq

    def cds_pos_to_omm_pos(self, cds_pos, ref_nuc):
        if self.error != "": return self.error

        cds_pos_ali = self.pos_add_gap(cds_pos, self.cds_seq_aligned)
        if self.cds_seq_aligned[cds_pos_ali] != ref_nuc:
            return "RefOMMDiff"

        if self.omm_seq_aligned[cds_pos_ali] == "-":
            return "NotConvertible"

        nbr_nuc_omm_ali = cds_pos_ali - self.omm_seq_aligned[:cds_pos_ali].count("-")
        omm_pos = self.pos_add_gap(nbr_nuc_omm_ali, self.omm_seq)

        if omm_pos is None or omm_pos >= len(self.omm_seq):
            return "NotConvertible"
        elif self.omm_seq[omm_pos] != ref_nuc:
            self.file_err.write(self.ensg + '\n')
            self.file_err.write("REF {0}\n".format(ref_nuc))
            self.file_err.write("CDS aligned position {0}:{1}\n".format(cds_pos_ali, self.cds_seq_aligned[cds_pos_ali]))
            self.file_err.write(str(self) + '\n')
            self.file_err.write("OMM aligned position {0}:{1}\n".format(cds_pos_ali, self.omm_seq_aligned[cds_pos_ali]))
            self.file_err.write("OMM position {0}:{1}\n".format(omm_pos, self.omm_seq[omm_pos]))
            return "RefOMMDiff"
        return omm_pos

    @staticmethod
    def score_ali(seq1, seq2):
        score = 0
        assert len(seq1) == len(seq2)
        for c1, c2 in zip(seq1, seq2):
            if c1 == c2:
                score += 1
        return score

    @staticmethod
    def pos_add_gap(pos_seq_ungap, seq_gap):
        tmp_pos_seq_ungap = -1
        for index, nuc in enumerate(seq_gap):
            if nuc != "-":
                tmp_pos_seq_ungap += 1
            if tmp_pos_seq_ungap == pos_seq_ungap:
                return index

    @staticmethod
    def chunks(l, n):
        """ Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    @staticmethod
    def gapsFromPeptide(peptide_seq, nucleotide_seq):
        if len(nucleotide_seq) == peptide_seq * 3: return nucleotide_seq
        codons = [codon for codon in Alignment.chunks(nucleotide_seq, 3)]  # splits nucleotides into codons (triplets)
        gapped_codons = []
        codon_count = 0
        for aa in peptide_seq:  # adds '---' gaps to nucleotide seq corresponding to peptide
            if aa != '-':
                gapped_codons.append(codons[codon_count])
                codon_count += 1
            else:
                gapped_codons.append('---')
        return ''.join(gapped_codons)


def extract_fasta(file_path, specie):
    if not os.path.exists(file_path):
        file_path = file_path.replace("_null_", "__")
    for f in SeqIO.parse(open(file_path, 'r'), 'fasta'):
        if f.id == specie:
            return str(f.seq)
    return ""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', required=True, type=str,
                        dest="f", metavar="<fasta>",
                        help="The relative name of the .fasta file")
    parser.add_argument('-g', '--gtf', required=True, type=str,
                        dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    parser.add_argument('-v', '--vcf', required=True, type=str,
                        dest="v", metavar="<vcf>",
                        help="The relative name of the .vcf file")
    parser.add_argument('-a', '--ali', required=True, type=str,
                        dest="ali", metavar="<ali>",
                        help="The alignment folder")
    parser.add_argument('-x', '--xml', required=True, type=str,
                        dest="xml", metavar="<xml>",
                        help="The xml folder")
    parser.add_argument('-s', '--species', required=True, type=str,
                        dest="species", metavar="<species>",
                        help="The species name")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest="output", metavar="<output>",
                        help="The output file")
    args = parser.parse_args()

    path = os.getcwd()

    dict_alignment = {}
    dict_cds, not_confirmed_tr = build_dict_cds(path, args.g)
    dict_tr_id = build_dict_trID(args.xml, args.species)

    print('Loading fasta file...')
    dict_fasta = {}
    for fasta in SeqIO.parse(gzip.open("{0}/{1}".format(path, args.f), 'rt'), 'fasta'):
        dict_fasta[fasta.id.split(".")[0]] = str(fasta.seq[:-3])
    print('Fasta file loaded.')

    annot_file = gzip.open(args.output, 'wt')
    error_file = open(args.output.replace(".vcf.gz", ".errors.tsv"), 'w')

    dict_cat_info = {"Syn": "{0} SNPs are synonymous variations",
                     "NonSyn": "{0} SNPs are non-synonymous variations",
                     "Stop": "{0} SNPs are stop variations",
                     "RefStop": "{0} SNPs have stop codon as reference amino-acid",
                     "RefDiff": "{0} SNPs retrieved from the fasta are not equal to the reference",
                     "NotIdentified": "{0} SNPs have non-identified reference or alternate amino-acid",
                     "NotInCds": "{0} SNPs are not inside the CDS",
                     "TrIdNotInOMM": "{0} transcript are not in OrthoMam",
                     "SeqOMMEmpty": "{0} fasta sequences from OrthoMam MSA are empy",
                     "SeqOMMnotMultiple3": "{0} fasta sequences from OrthoMam MSA are not multiple of 3",
                     "SeqCDSnotMultiple3": "{0} fasta sequences from CDS are not multiple of 3",
                     "NoAlignment": "{0} fasta sequences could not be aligned",
                     "NotConvertible": "{0} SNPs positions are not convertible",
                     "RefOMMDiff": "{0} SNPs retrieved from OMM are not equal to the reference"
                     }
    cat_errors = set(dict_cat_info.keys()).difference({"Syn", "NonSyn", "Stop"})
    dict_cat_nbr = {snp_cat: 0 for snp_cat in dict_cat_info}

    vcf_file = gzip.open("{0}/{1}".format(path, args.v), 'rt')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            suffix = "\tCHR\tTR_START\tTR_END\tTR_IDS\tENSG\tENSG_POS\tSNP_TYPE\tSNP_TYPES\n"
            annot_file.write(vcf_line.strip() + suffix)
            continue

        chromo, pos, snp_id, ref, alt = vcf_line.split("\t", maxsplit=5)[:5]
        if (ref not in nucleotides) or (alt not in nucleotides):
            continue

        snp_types = dict()
        transcript_id_list = vcf_line[vcf_line.rfind('\t') + 1:-1].split(",")
        for transcript_id in transcript_id_list:
            snp_type = dict_cds[transcript_id].snp_type(dict_fasta[transcript_id], int(pos), ref, alt)
            snp_types[transcript_id] = snp_type

        type_not_errors = {k: v for k, v in snp_types.items() if v not in cat_errors}
        if len(type_not_errors) == 0:
            dict_cat_nbr[most_common([v for v in snp_types.values() if v in cat_errors])] += 1
            continue

        ali_pos = dict()
        for tr_id in type_not_errors:
            if tr_id not in dict_tr_id:
                ali_pos[tr_id] = "TrIdNotInOMM"
                continue

            if tr_id not in dict_alignment:
                fasta_omm = extract_fasta("{0}{1}_NT.fasta".format(args.ali, dict_tr_id[tr_id]), args.species)
                dict_alignment[tr_id] = Alignment(dict_fasta[tr_id], fasta_omm, dict_tr_id[tr_id], error_file)

            nuc = complement[ref] if dict_cds[tr_id].strand == "-" else ref
            ali_pos[tr_id] = dict_alignment[tr_id].cds_pos_to_omm_pos(dict_cds[tr_id].nt_position(int(pos)), nuc)

        ali_pos_not_errors = {k: v for k, v in ali_pos.items() if v not in cat_errors}
        vcf_line = vcf_line.strip()
        if len(ali_pos_not_errors) == 0:
            dict_cat_nbr[most_common([v for v in ali_pos.values() if v in cat_errors])] += 1
            vcf_line += "\tNone\tNone"
        else:
            common_pos = most_common(list(ali_pos_not_errors.values()))
            vcf_line += "\t" + ",".join(set(ali_pos_not_errors.keys())) + "\t" + str(common_pos)

        types = [type_not_errors[k] for k in (ali_pos_not_errors if len(ali_pos_not_errors) > 0 else type_not_errors)]
        max_type = most_common(types)

        vcf_line += "\t" + max_type + "\t" + ",".join(types) + "\n"
        annot_file.write(vcf_line)
        dict_cat_nbr[max_type] += 1
    vcf_file.close()

    nbr_snp_total = sum(dict_cat_nbr.values())
    error_file.write("{0} SNPs in total".format(nbr_snp_total))
    for cat, nbr in dict_cat_nbr.items():
        error_file.write("\n\n" + dict_cat_info[cat].format(nbr) + " ({0:.3f}%)".format(nbr * 100. / nbr_snp_total))

    annot_file.close()
    error_file.close()

    print("{0} variants analyzed in total".format(nbr_snp_total))
    print("File containing {0} stop variants".format(dict_cat_nbr["Stop"]))
    print("File containing {0} synonymous variants".format(dict_cat_nbr["Syn"]))
    print("File containing {0} non-synonymous variants".format(dict_cat_nbr["NonSyn"]))
    nbr_errors = nbr_snp_total - dict_cat_nbr["Stop"] - dict_cat_nbr["Syn"] - dict_cat_nbr["NonSyn"]
    print("{0} errors variants".format(nbr_errors))
