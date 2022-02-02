import argparse
import os
from collections import defaultdict
from libraries_plot import *

# theta_label_dict = {"watterson": "$\\theta_W$ (Watterson)", "tajima": "$\\theta_{\\pi}$ (Tajima)",
#                    "fay_wu": "$\\theta_{H}$ (Fay and Hu)", "D_tajima": "Tajima's $d$", "H_fay_wu": "Fay and Hu's $H$"}
theta_label_dict = {"D_tajima": "Tajima's $D$", "H_fay_wu": "Fay and Hu's $H$"}

sorted_label = {"$1<S$": 6, "$0<S<1$": 5, "Synonymous": 4, "$-1<S<0$": 3, "$-3<S<-1$": 2, "$S<-3$": 1,
                "$0.8<SIFT$": 3, "$0.3<SIFT<0.8$": 2, "$0.1<SIFT<0.3$": 1, "$SIFT<0.1$": 0}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sfs', required=False, type=str, nargs="+", dest="tsv", help="Input pdf file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    parser.add_argument('-l', '--sample_list', required=False, type=str, dest="sample_list", help="Sample list file")

    args = parser.parse_args()
    concat_list = []
    for filepath in args.tsv:
        ddf = pd.read_csv(filepath.replace(".pdf", ".tsv"), sep="\t")
        species, pop, method = os.path.basename(filepath).replace(".pdf", "").split(".")
        ddf["pop"] = format_pop(pop.replace("_", " "))
        ddf["species"] = species
        ddf["method"] = method
        concat_list.append(ddf)

    df = pd.concat(concat_list)
    df["D_tajima"] = df["tajima"] - df["watterson"]
    df["H_fay_wu"] = df["tajima"] - df["fay_wu"]
    df.to_csv(args.output, sep="\t", index=False)

    for method, df_method in df.groupby(["method"]):
        for theta, theta_label in theta_label_dict.items():
            dico_theta = defaultdict(dict)
            pop2sp = {}
            for (pop, category), ddf in df_method.groupby(["pop", "category"]):
                assert len(ddf) == 1
                pop2sp[pop] = ddf["species"].values[0]
                dico_theta[pop][category] = ddf[theta].values[0]

            categories = set()
            for s in dico_theta.values():
                categories = categories.union(s.keys())
            categories = list(sorted(categories, key=lambda x: sorted_label[x]))

            species = [pop for pop, sp in sorted(pop2sp.items(), key=lambda kv: sp_sorted(*kv))]
            theta_matrix = np.ones((len(species), len(categories)))
            theta_matrix[:] = np.NaN

            for id_row, sp in enumerate(species):
                for id_col, category_set in enumerate(categories):
                    if sp in dico_theta and category_set in dico_theta[sp]:
                        theta_matrix[id_row, id_col] = dico_theta[sp][category_set]

            fig, ax = plt.subplots(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
            im, cbar = heatmap(theta_matrix.T, categories, species, ax=ax, cmap='RdBu_r', cbarlabel=theta_label)
            texts = annotate_heatmap(im, valfmt=lambda p: "{0:.2f}".format(p), div=True, fontsize=5)
            plt.tight_layout()
            plt.savefig(f"{args.output.replace('.tsv', '')}/{method}.{theta}.pdf", format="pdf")
            plt.clf()
            plt.close('all')
