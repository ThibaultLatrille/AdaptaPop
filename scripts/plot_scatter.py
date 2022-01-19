import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, Normalize
from scipy.ndimage.filters import gaussian_filter
import argparse
from libraries import *

RED = to_rgb("#EB6231")
BLUE = to_rgb("#5D80B4")
GREEN = to_rgb("#8FB03E")
GREY = to_rgb("grey")
BLACK = to_rgb("black")

error_kwargs = {"lw": .5, "zorder": -1}
my_dpi = 256
fontsize = 14
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)
idf = np.linspace(0, 1, 30)
plt.ylabel("$\\omega$", fontsize=fontsize)
plt.xlabel("$\\omega_0$", fontsize=fontsize)
xmin, xmax = 0.05, 1.0
ymin, ymax = 0.05, 1.0
plt.xlim((xmin, xmax))
plt.plot(idf, idf, color="black")


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return f"{s}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-g', '--granularity', required=True, type=str, dest="granularity", help="Gene or site level")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    gene = "gene" in args.granularity.lower()
    list_ensg = [i[:-3] for i in os.listdir(args.folder)]
    dico_omega_0, dico_omega = build_divergence_dico(args.folder, list_ensg, gene_level=gene)
    strg_ada_dico, ada_dico, nn_dico, unclass_dico = split_outliers(dico_omega_0, dico_omega, gene_level=gene)

    unclass_omega_0 = filtered_table_omega(dico_omega_0, unclass_dico, gene_level=gene)
    unclass_omega = filtered_table_omega(dico_omega, unclass_dico, gene_level=gene)
    strg_ada_omega_0 = filtered_table_omega(dico_omega_0, strg_ada_dico, gene_level=gene)
    strg_ada_omega = filtered_table_omega(dico_omega, strg_ada_dico, gene_level=gene)
    ada_omega_0 = filtered_table_omega(dico_omega_0, ada_dico, gene_level=gene)
    ada_omega = filtered_table_omega(dico_omega, ada_dico, gene_level=gene)
    nn_omega_0 = filtered_table_omega(dico_omega_0, nn_dico, gene_level=gene)
    nn_omega = filtered_table_omega(dico_omega, nn_dico, gene_level=gene)

    if gene:
        print('{0} adaptive genes'.format(len(ada_dico)))
        print('{0} nearly-neutral genes'.format(len(nn_dico)))
        print('{0} unclassified genes'.format(len(unclass_dico)))
        t = sum([len(d) for d in [strg_ada_dico, ada_dico, nn_dico, unclass_dico]])
        print(r'{0} total genes.'.format(t))
        # model = sm.OLS(list(nn_omega[:, 1]), sm.add_constant(list(nn_omega_0[:, 1])))
        # results = model.fit()
        # b, a = results.params[0:2]
        # plt.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
        #    float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"), color=GREEN)

        plt.errorbar(ada_omega_0[:, 1], ada_omega[:, 1],
                     xerr=[ada_omega_0[:, 1] - ada_omega_0[:, 0], ada_omega_0[:, 2] - ada_omega_0[:, 1]],
                     yerr=[ada_omega[:, 1] - ada_omega[:, 0], ada_omega[:, 2] - ada_omega[:, 1]],
                     alpha=0.4, label=r"${0}$ adaptive genes".format(len(ada_dico)), c=RED,
                     fmt='o', marker=None, mew=0, ecolor=RED, zorder=10, lw=.5, markersize=3.0)

        plt.errorbar(nn_omega_0[:, 1], nn_omega[:, 1],
                     xerr=[nn_omega_0[:, 1] - nn_omega_0[:, 0], nn_omega_0[:, 2] - nn_omega_0[:, 1]],
                     yerr=[nn_omega[:, 1] - nn_omega[:, 0], nn_omega[:, 2] - nn_omega[:, 1]],
                     alpha=0.4, label=r"${0}$ nearly-neutral genes".format(len(nn_dico)), c=GREEN,
                     fmt='o', marker=None, mew=0, ecolor=GREEN, zorder=5, lw=.5, markersize=3.0)

        plt.errorbar(unclass_omega_0[:, 1], unclass_omega[:, 1],
                     xerr=[unclass_omega_0[:, 1] - unclass_omega_0[:, 0],
                           unclass_omega_0[:, 2] - unclass_omega_0[:, 1]],
                     yerr=[unclass_omega[:, 1] - unclass_omega[:, 0], unclass_omega[:, 2] - unclass_omega[:, 1]],
                     alpha=0.4, label=r"${0}$ unclassified genes".format(len(unclass_dico)), c=GREY,
                     fmt='o', marker=None, mew=0, ecolor=GREY, zorder=0, lw=.5, markersize=3.0)

    else:
        xmin, xmax = 0.05, 1.4
        ymin, ymax = 0.05, 1.4
        xbins, ybins = 100, 100
        colors = np.zeros((ybins, xbins, 4))
        list_plot = list()
        list_plot.append((GREY, unclass_omega_0, unclass_omega))
        list_plot.append((RED, strg_ada_omega_0, strg_ada_omega))
        list_plot.append((GREEN, nn_omega_0, nn_omega))

        heatmaps = list()
        for COLOR, omega_0, omega in list_plot:
            heatmap, xedges, yedges = np.histogram2d(omega_0[:, 1], omega[:, 1], bins=(xbins, ybins),
                                                     range=[[xmin, xmax], [ymin, ymax]])
            extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

            heatmaps.append((COLOR, extent, heatmap))

        max_heatmap = max([np.max(h) for c, e, h in heatmaps])
        for COLOR, extent, heatmap in heatmaps:
            if COLOR == GREY:
                continue
            colors[..., 0:3] = COLOR[0:3]
            colors[..., -1] = Normalize(0.0, 1.5 * np.log(max_heatmap), clip=True)(np.log(heatmap.T))
            # mask = colors[..., -1] == 0.0
            # colors[..., -1] = colors[..., -1] * 0.9 + 0.1
            # colors[..., -1][mask] = 0.0
            plt.imshow(colors, extent=extent, origin="lower", aspect="auto", interpolation='nearest')
            data = gaussian_filter(heatmap.T, 1.25)
            set_levels = list(np.logspace(np.log10(25 if np.max(data) >= 25 else 5), np.log10(np.max(data)) - 0.15,
                                          3 if COLOR == BLUE else 5))
            cs = plt.contour(data, extent=extent, levels=set_levels, colors=[COLOR] * len(set_levels), linestyles='-')
            plt.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=10)

        plt.scatter(0, 0, label=r'{0} adaptive sites'.format(sum([len(v) for v in strg_ada_dico.values()])),
                    color=RED)
        plt.scatter(0, 0, label=r'{0} nearly-neutral sites'.format(sum([len(v) for v in nn_dico.values()])),
                    color=GREEN)
        plt.scatter(0, 0, label=r'{0} unclassifed sites'.format(sum([len(v) for v in unclass_dico.values()])),
                    color=GREY)
        print(r'{0} unclassified sites'.format(sum([len(v) for v in unclass_dico.values()])))
        t = sum([sum([len(v) for v in d.values()]) for d in [strg_ada_dico, ada_dico, nn_dico, unclass_dico]])
        print(r'{0} total sites.'.format(t))
    plt.legend(fontsize=fontsize_legend, loc="lower right")
    plt.ylim((ymin, ymax))
    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")

    print('Plot completed')
