from math import floor
import argparse
from Bio.Phylo.PAML import yn00
from Bio import SeqIO
from libraries import *

GREEN = "#8FB03E"
RED = "#EB6231"
dico_dNdS = {}


def dNdS_data_frame(ensg, specie_1, specie_2, ali_folder, fixed_poly, tmp_path):
    ali_seq = "{0}/{1}.fasta".format(tmp_path, ensg)
    seq_dico = {}
    for f in SeqIO.parse(open(ali_folder + "/" + ensg + ".fasta", 'r'), 'fasta'):
        if f.id in [specie_1, specie_2]:
            mutable_seq = f.seq.tomutable()
            for j in range(len(mutable_seq) // 3):
                codon = mutable_seq[j * 3:j * 3 + 3]
                if codontable[str(codon)] == "X":
                    mutable_seq[j * 3:j * 3 + 3] = "---"
            f.seq = mutable_seq
            seq_dico[f.id] = f

    if (specie_1 not in seq_dico) or (specie_2 not in seq_dico):
        dico_dNdS[ensg] = (0, 0, 0, 0)
        return

    if ensg in fixed_poly:
        for pos, (ref, alt) in fixed_poly[ensg].items():
            if seq_dico[specie_1].seq[pos] != ref:
                assert seq_dico[specie_1].seq[pos] == complement[ref]
                seq_dico[specie_1].seq[pos] = complement[alt]
            else:
                seq_dico[specie_1].seq[pos] = alt

    SeqIO.write([seq_dico[specie_1], seq_dico[specie_2]], ali_seq, "fasta")
    yn = yn00.Yn00(alignment=ali_seq, working_dir=tmp_path, out_file="{0}/{1}_yn.out".format(tmp_path, ensg))
    try:
        res = yn.run()[specie_1][specie_2]['YN00']
        dn = int(round(res['N'] * res['dN']))
        ds = int(round(res['S'] * res['dS']))
        dico_dNdS[ensg] = (res['N'], dn, res['S'], ds)
    except IndexError:
        dico_dNdS[ensg] = (0, 0, 0, 0)


def snp_data_frame(vcf_path):
    print("Loading file " + vcf_path)
    fixed_poly = dict()
    snp_table, header = {"ENSG": [], "POS": [], "TYPE": [], "COUNT": [], "SAMPLE_SIZE": []}, {}
    vcf_file = gzip.open(vcf_path, 'rt')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                header = {k: i for i, k in enumerate(vcf_line.strip().split("\t"))}
            continue

        line_list = vcf_line.strip().split("\t")
        ensg = line_list[header["ENSG"]]
        if ensg == "None": continue
        ensg += "_NT"

        if line_list[header["CHR"]] in ["X", "Y", "MT"]: continue

        genotypes = [s.split(":")[0].count("1") for s in line_list[header["FORMAT"] + 1:header["CHR"]] if
                     ("|" in s) or ("/" in s)]
        count = sum(genotypes)
        n = len(genotypes) * 2
        if count == 0: continue

        if count == n:
            if ensg not in fixed_poly: fixed_poly[ensg] = {}
            ref, alt = line_list[header["REF"]], line_list[header["ALT"]]
            fixed_poly[ensg][int(line_list[header["ENSG_POS"]])] = ref, alt
            continue

        snp_table["COUNT"].append(count)
        snp_table["ENSG"].append(ensg)
        snp_table["POS"].append(line_list[header["ENSG_POS"]])
        snp_table["TYPE"].append(line_list[header["SNP_TYPE"]])
        snp_table["SAMPLE_SIZE"].append(n)

    print(vcf_path + " loaded.")
    return pd.DataFrame(snp_table), fixed_poly


def dfe_alpha(filepath, df, ensg_list, species_focal, sister_species, ali_folder, fixed_poly, tmp_path):
    set_n = set(df["SAMPLE_SIZE"])
    assert len(set_n) == 1
    n = set_n.pop()

    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = np.zeros(n), np.zeros(n)
    for ensg in ensg_list:
        for row in df[df["ENSG"] == ensg].itertuples(index=False):
            if row.TYPE == "Syn":
                ps += 1
                sfs_s[row.COUNT] += 1
            elif row.TYPE == "NonSyn":
                pn += 1
                sfs_n[row.COUNT] += 1
            else:
                assert row.TYPE == "Stop"

        if ensg not in dico_dNdS:
            dNdS_data_frame(ensg, species_focal, sister_species, ali_folder, fixed_poly, tmp_path)

        sites_n += dico_dNdS[ensg][0]
        dn += dico_dNdS[ensg][1]
        sites_s += dico_dNdS[ensg][2]
        ds += dico_dNdS[ensg][3]

    sfs_list = ["Summed", n]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        range_sfs = range(1, int(floor(n // 2)) + 1)
        assert len(range_sfs) * 2 == n
        sfs_list += [nbr_site] + [(sfs[i] + sfs[n - i]) if n - i != i else sfs[i] for i in range_sfs]

    sfs_list += [sites_n, dn, sites_s, ds]
    content = "{0}+{1} ({2})".format(species_focal, sister_species, ",".join(ensg_list)) + "\n"
    content += "\t".join(map(str, sfs_list)) + "\n"
    if os.path.isfile(filepath):
        read = "".join(open(filepath, 'r').readlines())
        if read == content: return

    sfs_file = open(filepath, 'w')
    sfs_file.write(content)
    sfs_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--ali_folder', required=True, type=str, dest="ali_folder",
                        help="folder containing the OrthoMam alignments")
    parser.add_argument('-d', '--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-v', '--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('-g', '--gene_level', required=False, type=bool, default=True, dest="gene",
                        help="At the gene level")
    parser.add_argument('-f', '--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('-s', '--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-r', '--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('-t', '--tmp_folder', required=True, type=str, dest="tmp_folder",
                        help="Temporary files path")

    args = parser.parse_args()

    cds_list, table_omega_0, table_omega = build_divergence_tables(args.div_folder, gene_level=args.gene)
    cds_adaptive_list, cds_epistasis_list = detect_outliers(cds_list, table_omega_0, table_omega)
    cds_outliers_list = cds_adaptive_list + cds_epistasis_list
    non_outliers_list = [i for i in cds_list if cds_list not in cds_outliers_list]
    print('{0} outliers'.format(len(cds_outliers_list)))
    print('{0} non outliers'.format(len(cds_list) - len(cds_outliers_list)))

    df_snp, fix_poly = snp_data_frame(args.vcf)

    dfe_alpha("{0}/OUTLIERS.txt".format(args.output), df_snp, cds_adaptive_list, args.focal_species,
              args.sister_species, args.ali_folder, fix_poly, args.tmp_folder)

    np.random.seed(123456789)

    for rep in range(1, args.rep + 1):
        rand_set = np.random.choice(non_outliers_list, len(cds_adaptive_list), replace=False)

        dfe_alpha("{0}/{1}.txt".format(args.output, rep), df_snp, rand_set, args.focal_species,
                  args.sister_species, args.ali_folder, fix_poly, args.tmp_folder)
