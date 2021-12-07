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
    parser.add_argument('--ensembl', required=True, type=str, dest="ensembl", metavar="<ensembl>",
                        help="The path to the ensembl .vcf file")
    parser.add_argument('--nextgen', required=True, type=str, dest="nextgen", metavar="<nextgen>",
                        help="The path to the nextgen .vcf file")
    parser.add_argument('--dbSNP', required=True, type=str, dest="dbSNP", metavar="<dbSNP>",
                        help="The path to the dbSNP .vcf file")
    parser.add_argument('--output', required=True, type=str, dest="output", metavar="<output>",
                        help="The path to the output .vcf file")

    args = parser.parse_args()

    print("processing nextgen")
    header = ""
    dico_nextgen = defaultdict(dict)
    with gzip.open(args.nextgen, 'rt') as nextgen_file:
        for nextgen_line in nextgen_file:
            if nextgen_line[0] == "#":
                header += nextgen_line
                continue
            ng_sl = nextgen_line.split('\t')
            dico_nextgen[ng_sl[0]][int(ng_sl[1])] = (ng_sl[3], ng_sl[4], ng_sl[7:])
    print("nextgen processed")
    print(f"{sum([len(d) for d in dico_nextgen.values()])} SNPs are in NextGen.")

    print("processing dbSNP")
    dico_dbSNP = {}
    tot_dbSNP = 0
    with gzip.open(args.dbSNP, 'rt') as dbSNP_file:
        for dbSNP_line in dbSNP_file:
            if dbSNP_line[0] == "#":
                continue
            tot_dbSNP += 1
            if "SID=PRJEB6057" not in dbSNP_line:
                continue
            dbSNP_sl = dbSNP_line.split('\t')
            dico_dbSNP[dbSNP_sl[2]] = (dbSNP_sl[0].replace("chr", ""), int(dbSNP_sl[1]))
    print("dbSNP processed")
    dbSNP = len(dico_dbSNP)
    print(f"{dbSNP} SNPs out of {tot_dbSNP} from dbSNP are also in NextGen ({100 * dbSNP / tot_dbSNP:.2f}%).")

    print("processing ensembl")
    not_dbSNP, not_ng, not_ref, converted = 0, 0, 0, 0
    cpt, ref_chg, ref_chg_cpt, err = 0, 0, 0, 0
    outputfile = gzip.open(args.output, 'wt')
    outputfile.write(header)
    with gzip.open(args.ensembl, 'rt') as ensembl_file:
        for i, ensembl_line in enumerate(ensembl_file):
            if ensembl_line[0] == "#":
                continue
            ensembl_splitline = ensembl_line.split('\t')
            rsid = ensembl_splitline[2]
            if rsid not in dico_dbSNP:
                not_dbSNP += 1
                continue
            chromosome, pos = dico_dbSNP[rsid]
            if pos not in dico_nextgen[chromosome]:
                not_ng += 1
                continue
            ng_ref, ng_alt, info = dico_nextgen[chromosome][pos]
            if len(ng_ref) > 1 or len(ng_alt) > 1:
                not_ref += 1
                continue
            ref, alt = ensembl_splitline[3], ensembl_splitline[4]
            if ref != ng_ref:
                if ng_ref == complement[ref] and ng_alt == complement[alt]:
                    cpt += 1
                elif alt == ng_ref and ref == ng_alt:
                    ref_chg += 1
                    info = [convert_genotype(g) for g in info]
                elif alt == complement[ng_ref] and ref == complement[ng_alt]:
                    ref_chg_cpt += 1
                    info = [convert_genotype(g) for g in info]
                else:
                    err += 1
                    print("error")
                    print(rsid)
                    print(ng_ref, ng_alt, info)
                    print(ensembl_line)
                    continue
            outputline = "\t".join(ensembl_splitline[0:7]) + "\t" + "\t".join(info)
            outputfile.write(outputline)
            converted += 1

    print("ensembl processed")
    tot_ensembl = not_ng + not_dbSNP + not_ref + converted
    print(f"{not_dbSNP} SNPs out of {tot_ensembl} are not in dbSNP ({100 * not_dbSNP / tot_ensembl:.2f}%).")
    print(f"{not_ng} SNPs out of {tot_ensembl} are not in NextGen ({100 * not_ng / tot_ensembl:.2f}%).")
    print(f"{not_ref} SNPs out of {tot_ensembl} are not SNV ({100 * not_ref / tot_ensembl:.2f}%).")
    print(f"{err} SNPs out of {tot_ensembl} are errors ({100 * err / tot_ensembl:.2f}%).")
    print(f"{converted} successful SNPs liftover out of {tot_ensembl}  ({100 * converted / tot_ensembl:.2f}%).")

    print(f"{cpt} SNPs out of {converted} lifted changed cpt ({100 * cpt / converted:.2f}%).")
    print(f"{ref_chg} SNPs out of {converted} lifted changed ref ({100 * ref_chg / converted:.2f}%).")
    print(f"{ref_chg_cpt} SNPs out of {converted} lifted changed ref/cpt ({100 * ref_chg_cpt / converted:.2f}%).")
