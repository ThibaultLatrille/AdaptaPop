#!/usr/bin/env python3

import argparse
import os
from libraries import build_dict_cds

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gtf', required=True, type=str, dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", metavar="<output>",
                        help="The relative name of the output file")
    parser.add_argument('-c', '--chr2acc', required=False, type=str, default="", dest="chr2acc", metavar="<chr2acc>",
                        help="The relative name of the chr2acc file")
    parser.add_argument('-i', '--include', required=False, type=str, dest="i", metavar="<include>",
                        help="Information that must be contained in the CDS to be considered valid.")
    parser.add_argument('-p', '--chr_prefix', required=False, type=str, default="",
                        dest="chr_prefix", help="If the chromosome number is prefixed (e.g. chr1 instead of 1)")
    args = parser.parse_args()
    path = os.getcwd()

    dict_errors_info = {"NotConfirmed": "GTF file: {0} cds have a start or end that could not be confirmed",
                        "EmptyExon": "GTF file: {0} cds have at least one empty track",
                        "OverlappingExon": "GTF file: {0} cds have at least two exons overlapping",
                        "NotCoding": "GTF file: {0} cds are not a multiple of 3"}
    dict_errors_cds = {}
    for error in dict_errors_info:
        dict_errors_cds[error] = []

    dict_cds, not_confirmed_tr, pr2tr_id = build_dict_cds(path, args.g, args.chr2acc)

    bedfile = open("{0}/{1}".format(path, args.output), 'w')
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

        if cds.overlapping():
            dict_errors_cds["OverlappingExon"].append(transcript_id)
            continue

        for exon in cds.befile_lines(args.chr_prefix):
            bedfile.write(exon)

    bedfile.truncate()
    bedfile.close()

    nbr_cds_errors = sum([len(errors) for errors in dict_errors_cds.values()])
    nbr_cds_total = len(dict_cds)

    error_filename = '{0}.gtf_to_bed_errors.txt'.format(args.output.replace(".bed", ""))
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
