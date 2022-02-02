import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, Normalize

fontsize = 16
fontsize_legend = 14
my_dpi = 256
GREEN = "#8FB03E"
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
RED_RGB = to_rgb(RED)
BLUE_RGB = to_rgb(BLUE)
GREEN_RGB = to_rgb(GREEN)
GREY_RGB = to_rgb("grey")
BLACK_RGB = to_rgb("black")


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    im = ax.imshow(data, **kwargs)
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw, orientation='horizontal')
    cbar.ax.set_xlabel(cbarlabel)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
             rotation_mode="anchor")
    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    return im, cbar


def annotate_heatmap(im, data=None, div=False, valfmt="{x:.2f}", textcolors=("white", "black"), **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
    threshold_high = im.norm(data.max()) * (0.75 if div else 0.5)
    threshold_low = (im.norm(data.max()) * 0.25) if div else 0.0
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(threshold_high >= im.norm(data[i, j]) >= threshold_low)])
            text = im.axes.text(j, i, valfmt(data[i, j]), **kw)
            texts.append(text)
    return texts


def format_pval(p):
    return "0" if abs(p) < 1e-1 else "{0:.1g}".format(p)


def format_domega(p):
    return "{0:.2g}".format(p)


def column_format(data_frame):
    return "|" + "|".join(["l"] * 2 + ["r"] * (len(data_frame) - 2)) + "|"


def format_pop(t):
    if "up" == t:
        return "Equus"
    elif "dogs" == t:
        return "Canis"
    elif " " in t:
        return "".join([s[0] for s in t.split(" ")])
    else:
        return t


def sp_sorted(pop, sp):
    out = sp + "_" + pop
    if "Homo" in out:
        return "Z" + out
    elif "Chloro" in out:
        return "Y" + out
    elif "Ovis" in out:
        return "X" + out
    elif "Capra" in out:
        return "W" + out
    elif "Bos" in out:
        return "V" + out
    elif "Canis" in out:
        return "U" + out
    elif "Equus" in out:
        return "T" + out
    else:
        return out


def extend_pop(p, sample):
    fp = format_pop(p)
    if sample[p] == fp:
        return p
    else:
        return f"{sample[p]} ({fp})"


def sort_df(df, sample_list_path):
    df = df.iloc[df.apply(lambda r: sp_sorted(format_pop(r["pop"]), r["species"]), axis=1).argsort()]
    df["species"] = df.apply(lambda r: "Ovis orientalis" if r["pop"] == "IROO" else r["species"], axis=1)
    df["species"] = df.apply(lambda r: "Ovis vignei" if r["pop"] == "IROV" else r["species"], axis=1)
    df["species"] = df.apply(lambda r: "Capra aegagrus" if r["pop"] == "IRCA" else r["species"], axis=1)

    sample_iterrows = list(pd.read_csv(sample_list_path, sep='\t').iterrows())
    dico_sample = {r["SampleName"].replace("_", " "): r["Location"] for _, r in sample_iterrows}

    df["pop"] = df.apply(lambda r: extend_pop(r["pop"], dico_sample), axis=1)
    return df


def tex_f(x, highlight=False, pad=""):
    if x == 0:
        s = "0.0"
    elif 0.001 < abs(x) < 10:
        s = f"{x:6.3f}"
    elif 10 <= abs(x) < 10000:
        s = f"{x:6.1f}"
    else:
        s = f"{x:6.2g}"
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + '\\times 10^{' + str(int(exp)) + '}'
    if highlight:
        return "$\\bm{" + s + "{^*}}$"
    else:
        return f"${s}{pad}$"


def format_pval_df(d, prefix="", alpha=0.05):
    col = prefix + "pval_adj"
    d[col] = d.apply(lambda r: tex_f(r[col], r[col] < alpha, "~~"), axis=1)
    return d


def adjusted_holm_pval(d, prefix="", alpha=0.05, format_p=True):
    n = len(d[prefix + "pval"])
    sorted_pval = sorted(zip(d[prefix + "pval"], d.index))
    sorted_adjpval = [[min(1, pval * (n - i)), p] for i, (pval, p) in enumerate(sorted_pval)]
    for i in range(1, len(sorted_adjpval)):
        if sorted_adjpval[i][0] <= sorted_adjpval[i - 1][0]:
            sorted_adjpval[i][0] = sorted_adjpval[i - 1][0]
    holm = {p: pval for pval, p in sorted_adjpval}
    d[prefix + "pval_adj"] = [holm[p] for p in d.index]
    if format_p:
        d = format_pval_df(d, prefix=prefix, alpha=alpha)
    return d
