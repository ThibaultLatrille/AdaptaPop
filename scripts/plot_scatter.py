import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
from libraries import *

RED = "#EB6231"
BLUE = "#5D80B4"
GREEN = "#8FB03E"

my_dpi = 128
fontsize = 16
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)
idf = np.linspace(0, 1, 30)
plt.ylabel("$\\omega$", fontsize=fontsize)
plt.xlabel("$\\omega_0$", fontsize=fontsize)
xmin, xmax = 0.05, 1.0
ymin, ymax = 0.05, 1.0
plt.xlim((xmin, xmax))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-g', '--granularity', required=True, type=str, dest="granularity", help="Gene or site level")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, gene_level=gene)
    adaptive_dico, epistasis_dico, nearly_neutral_dico = split_outliers(dico_omega_0, dico_omega, gene_level=gene)

    if gene:
        print('{0} adaptive genes'.format(len(adaptive_dico)))
        print('{0} epistasis genes'.format(len(epistasis_dico)))
        print('{0} nearly-neutral genes'.format(len(nearly_neutral_dico)))
        nearly_neutral_omega_0 = filtered_table_omega(dico_omega_0, nearly_neutral_dico, gene_level=gene)
        nearly_neutral_omega = filtered_table_omega(dico_omega, nearly_neutral_dico, gene_level=gene)

        model = sm.OLS(nearly_neutral_omega, sm.add_constant(nearly_neutral_omega_0))
        results = model.fit()
        b, a = results.params[0:2]
        plt.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
            float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"), color=GREEN)

        plt.scatter(nearly_neutral_omega_0, nearly_neutral_omega, linewidth=3, edgecolors='none', alpha=0.85,
                    label=r"${0}$ nearly-neutral CDS".format(len(nearly_neutral_dico)), c=GREEN)

        plt.scatter(filtered_table_omega(dico_omega_0, adaptive_dico, gene_level=gene),
                    filtered_table_omega(dico_omega, adaptive_dico, gene_level=gene), edgecolors='none', alpha=0.85,
                    linewidth=3, label=r"${0}$ adaptive CDS".format(len(adaptive_dico)), c=RED)

        plt.scatter(filtered_table_omega(dico_omega_0, epistasis_dico, gene_level=gene),
                    filtered_table_omega(dico_omega, epistasis_dico, gene_level=gene), edgecolors='none', alpha=0.85,
                    linewidth=3, label=r"${0}$ epistasis CDS".format(len(epistasis_dico)), c=BLUE)
        plt.legend(fontsize=fontsize_legend)
    else:
        ymin, ymax = 0.05, 1.4
        plt.hist2d(table_omega(dico_omega_0, gene), table_omega(dico_omega, gene), bins=100,
                   range=[[xmin, xmax], [ymin, ymax]], norm=mpl.colors.LogNorm(),
                   cmap='Blues')

        plt.text(0.5, 0.5, '{0} nearly-neutral sites\n'.format(sum([len(v) for v in nearly_neutral_dico.values()])),
                 horizontalalignment='center', verticalalignment='center', fontweight="bold",
                 fontsize=fontsize_legend * 1.5, color=GREEN)
        plt.plot(idf, idf, color=GREEN)

        plt.text(0.1, 0.8, '{0} adaptive sites\n'.format(sum([len(v) for v in adaptive_dico.values()])),
                 horizontalalignment='left', verticalalignment='top', fontweight="bold",
                 fontsize=fontsize_legend * 1.5, color=RED)
        plt.plot(idf, idf + 0.1, color=RED)
        plt.plot(np.linspace(0, 0.9, 30), [1.0] * len(idf), color=RED)

        plt.text(0.8, 0.1, '{0} epistasis sites\n'.format(sum([len(v) for v in epistasis_dico.values()])),
                 horizontalalignment='right', verticalalignment='bottom', fontweight="bold",
                 fontsize=fontsize_legend * 1.5, color=BLUE)
        plt.plot(idf, idf - 0.1, color=BLUE)

    plt.ylim((ymin, ymax))
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")

    print('Plot completed')
