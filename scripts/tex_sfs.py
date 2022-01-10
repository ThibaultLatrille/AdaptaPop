import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--sfs', required=False, type=str, nargs="+", dest="sfs", help="Input sfs file (pdf)")
    parser.add_argument('-o', '--tex_doc', required=False, type=str, dest="tex_doc", help="Main document tex file")
    parser.add_argument('-t', '--tex_include', required=False, type=str, dest="tex_include", help="Include tex file")
    args = parser.parse_args()

    o = open(args.tex_include, 'w')
    for sfs in args.sfs:
        name = sfs.replace(".pdf", "").split("/")[-1].replace("_", " ").replace(".", " - ")
        o.write("\\subsection{" + name + "} \n \n")
        o.write("\\begin{minipage}{0.49\\linewidth} \n")
        o.write("\\includegraphics[width=\\linewidth, page=1]{" + sfs + "} \n")
        o.write("\\end{minipage}\n")
        o.write("\\begin{minipage}{0.49\\linewidth}\n")
        o.write("\\includegraphics[width=\\linewidth, page=1]{" + sfs.replace(".pdf", ".normalize.pdf") + "} \n")
        o.write("\\end{minipage}\n")
    o.close()

    tex_to_pdf = "pdflatex -synctex=1 -interaction=nonstopmode -output-directory={0} {1}".format(
        os.path.dirname(args.tex_doc), args.tex_doc)
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)
