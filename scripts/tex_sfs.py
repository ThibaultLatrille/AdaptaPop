import argparse
from glob import glob
from collections import defaultdict
import os
from libraries_plot import format_pop

dict_method = {"MutSel": "site-specific Mutation-Selection codon models.",
               "Omega": "site-specific codon models.",
               "Omega 0": "site-specific Mutation-Selection codon models (average over all amino acids).",
               "SIFT": "SIFT score",
               "WS": "weak to strong mutations (AT$\\rightarrow$GC)",
               "SW": "strong to weak mutations (GC$\\rightarrow$AT)"}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--sfs', required=False, type=str, nargs="+", dest="sfs", help="Input sfs file (pdf)")
    parser.add_argument('-o', '--tex_doc', required=False, type=str, dest="tex_doc", help="Main document tex file")
    parser.add_argument('-t', '--tex_include', required=False, type=str, dest="tex_include", help="Include tex file")
    args = parser.parse_args()

    heatmap_dict = defaultdict(list)
    for result in sorted(glob(f"{os.path.dirname(args.tex_doc)}/results/*.pdf")):
        method, theta = os.path.basename(result).replace(".pdf", "").replace("_", " ").split(".")
        heatmap_dict[method].append(result)

    o = open(args.tex_include, 'w')
    o.write("\\section{Tajima's D - Fay \\& Hu's H} \n \n")
    for method, list_heatmap in heatmap_dict.items():
        o.write("\\subsection{" + dict_method[method].capitalize() + '} \n')
        o.write("\\begin{center}\n")
        for heatmap in list_heatmap:
            o.write("\\includegraphics[width=0.75\\linewidth, page=1]{" + heatmap + "} \\\\\n")
        o.write("\\end{center}\n")

    nested_dict = defaultdict(lambda: defaultdict(dict))
    for sfs in sorted(args.sfs):
        sp, pop, method = os.path.basename(sfs).replace(".pdf", "").replace("_", " ").split(".")
        nested_dict[sp][pop][method] = sfs

    for sp, nested_dict_1 in nested_dict.items():
        o.write("\\section{" + sp + "} \n \n")
        for pop, nested_dict_2 in nested_dict_1.items():
            if " " in pop:
                o.write("\\subsection{" + f"{pop} ({format_pop(pop)})" + "} \n \n")
            else:
                o.write("\\subsection{" + format_pop(pop) + "} \n \n")

            for method, sfs in nested_dict_2.items():
                o.write("\\subsubsection*{SFS for " + dict_method[method] + '} \n')
                o.write("\\begin{minipage}{0.49\\linewidth} \n")
                o.write("\\includegraphics[width=\\linewidth, page=1]{" + sfs + "} \n")
                o.write("\\end{minipage}\n")
                o.write("\\begin{minipage}{0.49\\linewidth}\n")
                o.write(
                    "\\includegraphics[width=\\linewidth, page=1]{" + sfs.replace(".pdf", ".normalize.pdf") + "} \n")
                o.write("\\end{minipage}\n")
                o.write("\\\\ \n")
                if "SIFT" in method:
                    o.write("\\subsubsection*{SIFT score versus Mutation-Selection model} \n")
                    sfs = sfs.replace("SIFT", "SIFT_vs_MutSel")
                    o.write("\\includegraphics[width=0.5\\linewidth, page=1]{" + sfs + "} \n")
                    o.write("\\\\ \n")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {1}".format(
        os.path.dirname(args.tex_doc), args.tex_doc)
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
