from libraries import build_divergence_dico, load_alignments
import bz2
import _pickle as cPickle
import os
from lxml import etree
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--xml', required=True, type=str, dest="xml", help="The xml folder")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--ali_folder', required=True, type=str, dest="ali_folder",
                        help="folder containing the OrthoMam alignments")
    parser.add_argument('--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('--level', required=False, type=str, default="gene", dest="level",
                        help="Gene or site level")
    parser.add_argument('-p', '--pp', required=True, type=str, dest="pp", help="Posterior probability")
    parser.add_argument('--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    args = parser.parse_args()

    dico_trid = {"ENSG": [], "TRID": [], "CHR": [], "STRAND": [], }
    for file in os.listdir(args.xml):
        root = etree.parse(args.xml + "/" + file).getroot()
        dico_trid["ENSG"].append(file.replace(".xml", ""))
        trid, chr, strand = "None", "None", "None"
        for info in root.findall(".//infoCDS[@specy='{0}']".format(args.focal_species)):
            trid = str(info.find('ensidTr').text)
            chr = str(info.find('chr').text)
            strand = str(info.find('strand').text)

        dico_trid["TRID"].append(trid)
        dico_trid["CHR"].append(chr)
        dico_trid["STRAND"].append(strand)

    gene = args.level.lower() == "gene"
    ensg_annot = pd.DataFrame(dico_trid)
    filter_set = set(ensg_annot[-ensg_annot["CHR"].isin(["X", "Y", "MT", "None"])]["ENSG"])
    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, filter_set, gene_level=gene, pp=args.pp)
    dico_alignments = load_alignments(dico_omega, args.focal_species, args.sister_species, args.ali_folder)
    filter_set = filter_set.intersection(set(dico_alignments))

    with bz2.BZ2File(args.output, "w") as file:
        cPickle.dump((filter_set, dico_omega_0, dico_omega, dico_alignments), file)
