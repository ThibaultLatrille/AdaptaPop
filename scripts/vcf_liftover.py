#!/usr/bin/env python3
import gzip
import argparse
from libraries import complement
from collections import defaultdict


def convert_genotype(gen):
    if "0|0:" in gen:
        return gen.replace("0|0:", "1|1:")
    elif "0|1:" in gen:
        return gen.replace("0|1:", "1|0:")
    elif "1|0:" in gen:
        return gen.replace("1|0:", "0|1:")
    elif "1|1:" in gen:
        return gen.replace("1|1:", "0|0:")
    else:
        return gen


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--target', required=True, type=str, dest="target", metavar="<target>",
                        help="The path to the target .vcf file")
    parser.add_argument('--source', required=True, type=str, dest="source", metavar="<source>",
                        help="The path to the source .vcf file")
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", metavar="<vcf>",
                        help="The path to the .vcf file to liftover")
    parser.add_argument('--output', required=True, type=str, dest="output", metavar="<output>",
                        help="The path to the output .vcf file")

    args = parser.parse_args()

    print("processing source")
    dico_source = defaultdict(dict)
    with gzip.open(args.source, 'rt') as source_file:
        for source_line in source_file:
            if source_line[0] == "#":
                continue
            source_sl = source_line.split('\t')
            dico_source[source_sl[0].replace("chr", "")][int(source_sl[1])] = source_sl[2]
    print("source processed")
    source = sum([len(d) for d in dico_source.values()])
    print(f"{source} SNPs in target vcf.")

    print("processing target")
    dico_target = dict()
    tot_target = 0
    with gzip.open(args.target, 'rt') as target_file:
        for target_line in target_file:
            if target_line[0] == "#":
                continue
            tot_target += 1
            target_sl = target_line.split('\t')
            dico_target[target_sl[2]] = (target_sl[0], int(target_sl[1]), target_sl[3], target_sl[4])
    print("target processed")
    print(f"{len(dico_target)} SNPs in target vcf.")

    print("processing vcf")
    not_source, not_target, not_ref, converted = 0, 0, 0, 0
    cpt, ref_chg, ref_chg_cpt, err = 0, 0, 0, 0
    outputfile = gzip.open(args.output, 'wt')
    not_in_source_file = gzip.open(args.output.replace(".vcf.gz", ".Not_in_source.vcf.gz"), 'wt')
    not_in_target_file = gzip.open(args.output.replace(".vcf.gz", ".Not_in_target.vcf.gz"), 'wt')
    with gzip.open(args.vcf, 'rt') as vcf_file:
        for i, vcf_line in enumerate(vcf_file):
            if vcf_line[0] == "#":
                outputfile.write(vcf_line)
                not_in_source_file.write(vcf_line)
                not_in_target_file.write(vcf_line)
                continue
            vcf_sl = vcf_line.split('\t')
            vcf_chr, vcf_pos, vcf_ref, vcf_alt = vcf_sl[0], int(vcf_sl[1]), vcf_sl[3], vcf_sl[4]
            if len(vcf_ref) > 1 or len(vcf_alt) > 1:
                not_ref += 1
                continue
            if vcf_pos not in dico_source[vcf_chr]:
                not_in_source_file.write(vcf_line)
                not_source += 1
                continue
            rsid = dico_source[vcf_chr][vcf_pos]
            if rsid not in dico_target:
                not_in_target_file.write(vcf_line)
                not_target += 1
                continue
            target_chr, target_pos, ref, alt = dico_target[rsid]
            info = vcf_sl[7:]
            if ref != vcf_ref:
                if vcf_ref == complement[ref] and vcf_alt == complement[alt]:
                    cpt += 1
                elif alt == vcf_ref and ref == vcf_alt:
                    ref_chg += 1
                    info = [convert_genotype(g) for g in info]
                elif alt == complement[vcf_ref] and ref == complement[vcf_alt]:
                    ref_chg_cpt += 1
                    info = [convert_genotype(g) for g in info]
                else:
                    err += 1
                    print("error")
                    print(rsid)
                    print(vcf_ref, vcf_alt, info)
                    print(vcf_line)
                    continue
            vcf_sl[0], vcf_sl[1], vcf_sl[3], vcf_sl[4] = target_chr, str(target_pos), ref, alt
            outputline = "\t".join(vcf_sl[0:7]) + "\t" + "\t".join(info)
            outputfile.write(outputline)
            converted += 1
    outputfile.close()
    not_in_source_file.close()
    not_in_target_file.close()
    print("target processed")
    tot_vcf = not_target + not_source + not_ref + converted
    print(f"{not_source} SNPs out of {tot_vcf} are not in source ({100 * not_source / tot_vcf:.2f}%).")
    print(f"{not_target} SNPs out of {tot_vcf} are not in target ({100 * not_target / tot_vcf:.2f}%).")
    print(f"{not_ref} SNPs out of {tot_vcf} are not SNV ({100 * not_ref / tot_vcf:.2f}%).")
    print(f"{err} SNPs out of {tot_vcf} are errors ({100 * err / tot_vcf:.2f}%).")
    print(f"{converted} successful SNPs liftover out of {tot_vcf}  ({100 * converted / tot_vcf:.2f}%).")

    print(f"{cpt} SNPs out of {converted} lifted changed cpt ({100 * cpt / converted:.2f}%).")
    print(f"{ref_chg} SNPs out of {converted} lifted changed ref ({100 * ref_chg / converted:.2f}%).")
    print(f"{ref_chg_cpt} SNPs out of {converted} lifted changed ref/cpt ({100 * ref_chg_cpt / converted:.2f}%).")
