import argparse
import os
from libraries import build_divergence_dico, split_outliers
from libraries_plot import plt
from upsetplot import UpSet, from_contents
from venn import venn

fontsize = 16
fontsize_legend = 14
my_dpi = 256
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--granularity', required=False, type=str, nargs="+", dest="granularity")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output folder")
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    args = parser.parse_args()

    list_ensg = [i[:-3] for i in os.listdir(args.folder)]

    dico_set, dico_div = {}, {}
    for granularity in args.granularity:
        level, method, pp = granularity.split("-")
        gene = "gene" in level

        if (gene, pp) not in dico_div:
            dico_div[(gene, pp)] = build_divergence_dico(args.folder, list_ensg, gene_level=gene, pp=pp)

        omega_0, omega = dico_div[(gene, pp)]
        adap_dico, _, _ = split_outliers(omega_0, omega, gene_level=gene, method=method)
        if gene:
            cds_focal_set = set(adap_dico)
        else:
            cds_focal_set = set()
            for ensg, site_list in adap_dico.items():
                if len(site_list) > 1:
                    cds_focal_set.add(ensg)

        dico_set[granularity.replace("_", " ")] = cds_focal_set

    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    plt.subplot(1, 1, 1)
    upset_dataframe = from_contents(dico_set)
    UpSet(upset_dataframe, subset_size='count', show_counts=True, min_subset_size=0, sort_by="cardinality").plot()
    plt.savefig(args.output.replace("venn_diagram", "upsetplot"), format="pdf")
    plt.clf()
    plt.close("all")

    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    plt.subplot(1, 1, 1)
    venn(dico_set)
    plt.savefig(args.output, format="pdf")
