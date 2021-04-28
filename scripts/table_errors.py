import argparse
import os
import pandas as pd


def find_errors(folder, gene_level=True, ci="0.025"):
    print('Loading divergence results.')
    ci = "0.0025" if gene_level else "0.05"
    list_errors = []
    for file in sorted(os.listdir(folder)):
        sitemutsel_path = "{0}/{1}/sitemutsel_1.run.ci{2}.tsv".format(folder, file, ci)
        siteomega_path = "{0}/{1}/siteomega_1.run.ci{2}.tsv".format(folder, file, ci)
        if not os.path.isfile(siteomega_path) or not os.path.isfile(sitemutsel_path):
            list_errors.append(file)

    return list_errors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--div_folder', required=False, default="OrthoMam", type=str, dest="div_folder",
                        help="folder containing OrthoMam results")

    args = parser.parse_args()
    df = pd.DataFrame({"ENSG": find_errors(args.div_folder, gene_level=True)})
    df.to_csv("TableErrors.tsv", sep="\t", index=False, header=False)
