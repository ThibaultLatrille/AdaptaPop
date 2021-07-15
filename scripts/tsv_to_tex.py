import argparse
import pandas as pd
import os
from libraries import tex_f

header = ["Species", "Population", "SFS", "Level", "Model",
          "$\\omega_{\\textrm{A}}^{\\textrm{S}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{S}}$",
          "$\\omega^{\\textrm{S}}$", "$\\alpha^{\\textrm{S}}$",
          "$\\omega_{\\textrm{A}}^{\\textrm{N}}$", "$\\omega_{\\textrm{NA}}^{\\textrm{N}}$",
          "$\\omega^{\\textrm{N}}$", "$\\alpha^{\\textrm{N}}$",
          "p-value", "$a$", "$r^2$"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-o', '--output', required=False, type=str, dest="output", help="Output tex file")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df = df.groupby(["species", "pop", "sfs", "granularity", "model"]).max().reset_index()
    df.sort_values(by=["species", "pop", "granularity", "sfs", "model"], inplace=True)

    o = open(args.output, 'w')
    o.write(df.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=header))
    o.write("\\newpage")

    for model in ["grapes", "polyDFE", "dfem"]:
        ddf = df[df["model"] == model].drop(["model"], axis=1)
        if len(ddf["species"]) == 0: continue
        o.write("\\section*{" + model + "} \n")
        o.write(ddf.to_latex(index=False, escape=False, longtable=True, float_format=tex_f,
                             header=[i for i in header if i != "Model"]))
        o.write("\\newpage")

    sub_header = [i for i in header if i not in ["SFS", "Level", "Model"]]
    for sfs in ["folded", "unfolded"]:
        for granularity in ["gene", "site"]:
            for model in ["grapes", "polyDFE", "dfem"]:
                ddf = df[(df["sfs"] == sfs) & (df["granularity"] == granularity) & (df["model"] == model)].drop(
                    ["sfs", "granularity", "model"], axis=1)
                if len(ddf["species"]) == 0: continue
                o.write("\\section*{" + "{0} SFS at {1} level - {2}".format(sfs.upper(), granularity, model) + "} \n")
                o.write(ddf.to_latex(index=False, escape=False, longtable=True, float_format=tex_f, header=sub_header))
                o.write("\\newpage")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {0}/main-table.tex".format(
        os.path.dirname(args.output))
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
