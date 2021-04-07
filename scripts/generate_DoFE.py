from math import floor
import argparse
from Bio.Phylo.PAML import yn00
from Bio import SeqIO, Seq
from libraries import *
from scipy.stats import norm

GREEN = "#8FB03E"
RED = "#EB6231"


def snp_data_frame(vcf_path):
    print("Loading file " + vcf_path)
    fixed_poly = dict()
    sample_size_set = set()
    snp_table, header = {"ENSG": [], "POS": [], "TYPE": [], "COUNT": []}, {}
    vcf_file = gzip.open(vcf_path, 'rt')
    for vcf_line in vcf_file:
        if vcf_line[0] == '#':
            if vcf_line[1] != '#':
                header = {k: i for i, k in enumerate(vcf_line.strip().split("\t"))}
            continue

        line_list = vcf_line.strip().split("\t")
        ensg = line_list[header["ENSG"]]
        if ensg == "None": continue

        if line_list[header["CHR"]] in ["X", "Y", "MT"]: continue

        genotypes = [s.split(":")[0].count("1") for s in line_list[header["FORMAT"] + 1:header["CHR"]] if
                     ("|" in s) or ("/" in s)]
        count = sum(genotypes)
        n = len(genotypes) * 2
        if count == 0: continue

        nuc_pos = int(line_list[header["ENSG_POS"]])
        codon_pos = int(nuc_pos / 3)
        if count == n:
            if ensg not in fixed_poly: fixed_poly[ensg] = {}
            ref, alt = line_list[header["REF"]], line_list[header["ALT"]]
            fixed_poly[ensg][codon_pos] = nuc_pos % 3, ref, alt
            continue

        snp_table["COUNT"].append(count)
        snp_table["ENSG"].append(ensg)
        snp_table["POS"].append(codon_pos)
        snp_table["TYPE"].append(line_list[header["SNP_TYPE"]])
        sample_size_set.add(n)

    print(vcf_path + " loaded.")
    assert len(sample_size_set) == 1
    n = sample_size_set.pop()
    return pd.DataFrame(snp_table).groupby("ENSG"), fixed_poly, n


def load_alignments(ensg_list, sp_1, sp_2, ali_folder):
    ali_dico = dict()

    for ensg in ensg_list:
        sister_focal_seqs = dict()
        for f in SeqIO.parse(open(ali_folder + "/" + ensg + "_NT.fasta", 'r'), 'fasta'):
            if f.id not in [sp_1, sp_2]: continue
            sister_focal_seqs[f.id] = f.seq

        if sp_1 in sister_focal_seqs and sp_2 in sister_focal_seqs:
            ali_dico[ensg] = sister_focal_seqs

    return ali_dico


def filter_positions(sister_focal_seqs, sp_focal, fixed_poly, gene_level, positions):
    seq_dico = dict()

    for sp_id, seq in sister_focal_seqs.items():

        s_list = list()
        for pos in (range(len(seq) // 3) if gene_level else positions):
            codon = str(seq[pos * 3:pos * 3 + 3])

            if pos in fixed_poly and sp_id == sp_focal:
                frame, ref, alt = fixed_poly[pos]
                if codon[frame] != ref:
                    assert codon[frame] == complement[ref]
                    alt = complement[alt]

                codon = codon[:frame] + alt + codon[frame + 1:]

            if codontable[codon] == "X": codon = "---"
            s_list.append(codon)

        seq_dico[sp_id] = "".join(s_list)

    return seq_dico


def run_yn00(seq1, seq2, tmp_path, filepath):
    SeqIO.write([seq1, seq2], filepath + ".fasta", "fasta")
    yn = yn00.Yn00(alignment=filepath + ".fasta", working_dir=tmp_path, out_file=filepath + "_yn00.out")

    try:
        res = yn.run()[seq1.id][seq2.id]['YN00']
        dn = int(round(res['N'] * res['dN']))
        ds = int(round(res['S'] * res['dS']))
        return res['N'], dn, res['S'], ds
    except:
        return 0, 0, 0, 0


def dfe_alpha(filepath, df, n, ensg_dico_pos, gene_level, sp_1, sp_2, ali_dico, fixed_poly, tmp_path):
    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    s1, s2 = [], []

    for ensg in ensg_dico_pos:
        seq_dico = filter_positions(ali_dico[ensg], sp_1, fixed_poly[ensg] if ensg in fixed_poly else {},
                                    gene_level, ensg_dico_pos[ensg])

        s1.append(seq_dico[sp_1])
        s2.append(seq_dico[sp_2])

        if ensg not in df.groups: continue

        dff = df.get_group(ensg)
        if not gene_level:
            filt = dff["POS"].isin(ensg_dico_pos[ensg])
            dff = dff[filt]

        for row in dff.itertuples(index=False):
            if row.TYPE == "Syn":
                ps += 1
                sfs_s[row.COUNT] += 1
            elif row.TYPE == "NonSyn":
                pn += 1
                sfs_n[row.COUNT] += 1
            else:
                assert row.TYPE == "Stop"

    seq1 = SeqIO.SeqRecord(id=sp_1, name="", description="", seq=Seq.Seq("".join(s1)))
    seq2 = SeqIO.SeqRecord(id=sp_2, name="", description="", seq=Seq.Seq("".join(s2)))
    sites_n, dn, sites_s, ds = run_yn00(seq1, seq2, tmp_path, filepath.replace(".txt", ""))

    sfs_list = ["Summed", n]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        range_sfs = range(1, int(floor(n // 2)) + 1)
        assert len(range_sfs) * 2 == n
        sfs_list += [nbr_site] + [(sfs[i] + sfs[n - i]) if n - i != i else sfs[i] for i in range_sfs]

    sfs_list += [sites_n, dn, sites_s, ds]
    content = "{0}+{1} ({2} sites)".format(sp_1, sp_2, len(seq1.seq)) + "\n"
    content += "\t".join(map(str, sfs_list)) + "\n"
    if os.path.isfile(filepath):
        read = "".join(open(filepath, 'r').readlines())
        if read == content: return

    sfs_file = open(filepath, 'w')
    sfs_file.write(content)
    sfs_file.close()


def subsample_sites(cds_dico, nbr_sites, weights):
    dico_list = list(cds_dico)
    rand_dico = dict()
    interval = [0]
    for k, v in cds_dico.items():
        interval.append(interval[-1] + len(v))

    assert nbr_sites <= interval[-1]
    site_choices = sorted(np.random.choice(interval[-1], nbr_sites, replace=False, p=weights))

    insert = np.searchsorted(interval, site_choices, side='right') - 1

    for ins, site in zip(insert, site_choices):
        ensg = dico_list[ins]
        pos = cds_dico[ensg][site - interval[ins]]
        if ensg not in rand_dico: rand_dico[ensg] = []
        rand_dico[ensg].append(pos)
    return rand_dico


def subsample_genes(cds_dico, nbr_sites, weights):
    return {k: None for k in np.random.choice(list(cds_dico), nbr_sites, replace=False, p=weights)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--ali_folder', required=True, type=str, dest="ali_folder",
                        help="folder containing the OrthoMam alignments")
    parser.add_argument('-d', '--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('-v', '--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('-g', '--granularity', required=False, type=str, default=True, dest="granularity",
                        help="Gene or site level")
    parser.add_argument('-n', '--nbr_sites', required=False, type=int, default=50000, dest="nbr_sites",
                        help="Number of sites to subsample")
    parser.add_argument('-c', '--nbr_genes', required=False, type=int, default=150, dest="nbr_genes",
                        help="Number of genes to subsample")
    parser.add_argument('-f', '--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('-s', '--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('-w', '--weighted', required=False, type=str, default='False', dest="weighted",
                        help="Weighted sampling such as to control for omega")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-r', '--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('-t', '--tmp_folder', required=True, type=str, dest="tmp_folder",
                        help="Temporary files path")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, gene_level=gene)
    dico_alignments = load_alignments(dico_omega, args.focal_species, args.sister_species, args.ali_folder)
    adaptive_dico, epistasis_dico, nearly_neutral_dico = split_outliers(dico_omega_0, dico_omega, gene_level=gene,
                                                                        filter_set=set(dico_alignments))
    if gene:
        if args.nbr_genes == -1 or args.nbr_genes > len(adaptive_dico): args.nbr_genes = len(adaptive_dico)
        print('{0} adaptive genes'.format(len(adaptive_dico)))
        print('{0} epistasis genes'.format(len(epistasis_dico)))
        print('{0} nearly-neutral genes'.format(len(nearly_neutral_dico)))
    else:
        adapta_sites = sum([len(v) for v in adaptive_dico.values()])
        if args.nbr_sites == -1 or args.nbr_sites > adapta_sites: args.nbr_sites = adapta_sites
        print('{0} adaptive sites'.format(adapta_sites))
        print('{0} epistasis sites'.format(sum([len(v) for v in epistasis_dico.values()])))
        print('{0} nearly-neutral sites'.format(sum([len(v) for v in nearly_neutral_dico.values()])))

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf)
    np.random.seed(123456789)

    weigths_nearly_neutral = None
    if args.weighted.lower() == "true":
        adaptive_table = filtered_table_omega(dico_omega, adaptive_dico, gene)
        weigths_nearly_neutral = norm.pdf(filtered_table_omega(dico_omega, nearly_neutral_dico, gene),
                                          loc=np.mean(adaptive_table), scale=np.std(adaptive_table))
        weigths_nearly_neutral /= np.sum(weigths_nearly_neutral)

    for rep in range(1, args.rep + 1):
        if gene:
            adaptive_set = subsample_genes(adaptive_dico, args.nbr_genes, None)
            nearly_neutral_set = subsample_genes(nearly_neutral_dico, args.nbr_genes, weigths_nearly_neutral)
        else:
            adaptive_set = subsample_sites(adaptive_dico, args.nbr_sites, None)
            nearly_neutral_set = subsample_sites(nearly_neutral_dico, args.nbr_sites, weigths_nearly_neutral)

        dfe_alpha("{0}/ADAPTIVE_{1}.txt".format(args.output, rep), df_snp, sample_size, adaptive_set, gene,
                  args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder)

        dfe_alpha("{0}/NEARLY_NEUTRAL_{1}.txt".format(args.output, rep), df_snp, sample_size, nearly_neutral_set, gene,
                  args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder)
