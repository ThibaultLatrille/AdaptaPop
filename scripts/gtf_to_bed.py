#!/usr/bin/env python3

import argparse
import os
import gzip


class Cds(object):
    def __init__(self, chromosome, strand, name):
        self.chromosome = chromosome
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


def build_dict_cds(data_path, file_name, includes):
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gtf', required=True, type=str,
                        dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    parser.add_argument('-i', '--include', required=False, type=str,
                        dest="i", metavar="<include>",
                        help="Information that must be contained in the CDS to be considered valid.")
    args = parser.parse_args()
    path = os.getcwd()

    dict_errors_info = {"NotConfirmed": "GTF file: {0} cds have a start or end that could not be confirmed",
                        "EmptyExon": "GTF file: {0} cds have at least one empty track",
                        "NotCoding": "GTF file: {0} cds are not a multiple of 3"}
    dict_errors_cds = {}
    for error in dict_errors_info:
        dict_errors_cds[error] = []

    dict_cds, not_confirmed_tr = build_dict_cds(path, args.g, args.i)

    bedfile = open("{0}/{1}.bed".format(path, args.g.replace(".gtf.gz", "")), 'w')
    for transcript_id, cds in dict_cds.items():

        if transcript_id in not_confirmed_tr:
            dict_errors_cds["NotConfirmed"].append(transcript_id)
            continue

        if cds.empty_exon():
            dict_errors_cds["EmptyExon"].append(transcript_id)
            continue

        if cds.not_coding():
            dict_errors_cds["NotCoding"].append(transcript_id)
            continue

        for exon in cds.befile_lines():
            bedfile.write(exon)

    bedfile.truncate()
    bedfile.close()

    nbr_cds_errors = sum([len(errors) for errors in dict_errors_cds.values()])
    nbr_cds_total = len(dict_cds)

    error_filename = '{0}.gtf_to_bed_errors.txt'.format(args.g.replace(".gtf.gz", ""))
    error_header = "{0} errors out of {1} coding sequences ({2:.3f}%)".format(
        nbr_cds_errors, nbr_cds_total, nbr_cds_errors * 100. / nbr_cds_total)
    print(error_header)
    print("Errors written in {0}".format(error_filename))

    error_file = open('{0}/{1}'.format(path, error_filename), 'w')
    error_file.write(error_header)
    if nbr_cds_errors > 0:
        for error, list_cds in dict_errors_cds.items():
            nbr_cds = len(list_cds)
            if nbr_cds > 0:
                error_file.write("\n\n" + dict_errors_info[error].format(nbr_cds))
                error_file.write("\n" + "\t".join(list_cds))
    error_file.close()
