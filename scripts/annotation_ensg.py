import scipy.stats as st
from libraries import build_divergence_dico, split_outliers, tex_f
import os
from lxml import etree
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-x', '--xml', required=True, type=str, dest="xml", help="The xml folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-s', '--species', required=True, type=str, dest="species", help="The species name")
    args = parser.parse_args()

    dico_trid = {"ENSG": [], "TRID": [], "CHR": [], "STRAND": [],}
    for file in os.listdir(args.xml):
        root = etree.parse(args.xml + "/" + file).getroot()
        dico_trid["ENSG"].append(file.replace(".xml", ""))
        trid, chr, strand = "None", "None", "None"
        for info in root.findall(".//infoCDS[@specy='{0}']".format(args.species)):
            trid = str(info.find('ensidTr').text)
            chr = str(info.find('chr').text)
            strand = str(info.find('strand').text)

        dico_trid["TRID"].append(trid)
        dico_trid["CHR"].append(chr)
        dico_trid["STRAND"].append(strand)
    pd.DataFrame(dico_trid).to_csv(args.output, sep="\t", index=False)
    print('TR_ID to ENSG conversion done.')
