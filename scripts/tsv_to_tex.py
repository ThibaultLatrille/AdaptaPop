import argparse
import pandas as pd
import os


def format_float(x):
    if 0.001 < abs(x) < 10:
        return "{:6.3f}".format(x)
    elif 10 <= abs(x) < 10000:
        return "{:6.1f}".format(x)
    else:
        s = "{:6.2g}".format(x)
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + '$\\times 10^{' + str(int(exp)) + '}$'
            s = " " * (5 + 6 - len(s)) + s
        return s


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df = df.groupby(["species", "pop", "sfs", "granularity"]).max().reset_index()
    df.sort_values(by=["species", "pop", "granularity", "sfs"], inplace=True)
    df.to_latex(args.output, index=False, escape=False, longtable=True, float_format=format_float,
                header=["Species", "Population", "SFS", "Level", "$\\omega_{\\textrm{A}}^{\\textrm{adaptive}}$",
                        "$\\omega_{\\textrm{A}}^{\\textrm{nearly-neutral}}$", "p-value", "$a$", "$r^2$"])
    os.system("pdflatex -output-directory={0} {0}/main-table.tex".format(os.path.dirname(args.output)))
