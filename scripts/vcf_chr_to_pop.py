#!/usr/bin/env python3
import gzip
import pandas as pd
import argparse


def to_int(x):
    try:
        return int(x)
    except:
        return 99


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, nargs="+",
                        dest="vcf", metavar="<vcf>",
                        help="The path to the .vcf files")
    parser.add_argument('--panel', required=True, type=str,
                        dest="panel", metavar="<panel>",
                        help="The relative name of the .panel file")
    parser.add_argument('--output', required=True, type=str,
                        dest="output", metavar="<output>",
                        help="The path to the output .vcf files")

    args = parser.parse_args()
    panel = pd.read_csv(args.panel, sep='\t', usecols=('sample', 'pop'))
    individuals = set(panel["sample"])

    sorted_chr = sorted(args.vcf, key=lambda x: to_int(x.split(".")[-3].replace("chr", "")))

    for pop in sorted(set(panel['pop'])):
        pop_individuals = set(panel[panel["pop"] == pop]["sample"])
        print("{0} population with {1} individuals".format(pop, len(pop_individuals)))

        vcf_pop_file = gzip.open(args.output.replace("POPULATION", pop), 'wt')

        header_written = False
        for vcf_chr_path in sorted_chr:
            if 'wgs' in vcf_chr_path: continue
            print("Chromosome " + vcf_chr_path)
            vcf_chr_file = gzip.open(vcf_chr_path, 'rt')
            for line in vcf_chr_file:
                if not line.startswith('##'): break

            assert line.startswith('#')

            if not header_written:
                header = line.strip().split('\t')
                samples = pop_individuals.intersection(header)
                include_index = [cell in samples or cell not in individuals for cell in header]

                vcf_pop_file.write("\t".join([cell for i, cell in enumerate(header) if include_index[i]]) + "\n")
                header_written = True
            else:
                assert line.strip().split('\t') == header

            for line in vcf_chr_file:

                cells = []
                for i, cell in enumerate(line.strip().split('\t')):
                    if i >= len(header) or include_index[i]:
                        cells.append(cell)

                vcf_pop_file.write("\t".join(cells) + "\n")

            vcf_chr_file.close()

        vcf_pop_file.close()

    print("Analysis completed")
