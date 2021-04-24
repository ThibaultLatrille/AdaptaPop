import argparse
from libraries import *
from scipy.stats import norm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ali_folder', required=True, type=str, dest="ali_folder",
                        help="folder containing the OrthoMam alignments")
    parser.add_argument('--div_folder', required=True, type=str, dest="div_folder",
                        help="folder containing OrthoMam results")
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--granularity', required=False, type=str, default=True, dest="granularity",
                        help="Gene or site level")
    parser.add_argument('--nbr_sites', required=False, type=int, default=50000, dest="nbr_sites",
                        help="Number of sites to subsample")
    parser.add_argument('--nbr_genes', required=False, type=int, default=150, dest="nbr_genes",
                        help="Number of genes to subsample")
    parser.add_argument('--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('--weighted', required=False, type=str, default='False', dest="weighted",
                        help="Weighted sampling such as to control for omega")
    parser.add_argument('--seed', required=False, type=int, default=123456789, dest="seed",
                        help="Seed for random generator")
    parser.add_argument('--dfem_path', required=True, type=str, dest="dfem_path", help="Executable path for DFEM")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, gene_level=gene)
    dico_alignments = load_alignments(dico_omega, args.focal_species, args.sister_species, args.ali_folder)
    _, adaptive_dico, epistasis_dico, nearly_neutral_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=gene,
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
    np.random.seed(args.seed)

    weigths_nearly_neutral = None
    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for {3} {4} set\n"
    f = open(args.output + ".txt", 'w')
    w = filtered_table_omega(dico_omega, nearly_neutral_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, nearly_neutral_dico, gene)[:, 1]
    w_A = w - w_0
    f.write(txt.format("", np.mean(w), np.std(w), "all nearly-neutral", args.granularity))
    f.write(txt.format("0", np.mean(w_0), np.std(w_0), "all nearly-neutral", args.granularity))
    f.write(txt.format("A", np.mean(w_A), np.std(w_A), "all nearly-neutral", args.granularity))

    w = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, adaptive_dico, gene)[:, 1]
    w_A = w - w_0
    f.write(txt.format("", np.mean(w), np.std(w), "all adaptive", args.granularity))
    f.write(txt.format("0", np.mean(w_0), np.std(w_0), "all adaptive", args.granularity))
    f.write(txt.format("A", np.mean(w_A), np.std(w_A), "all adaptive", args.granularity))

    if args.weighted.lower() == "true":
        print("Weighted sampling controlling for w")
        adaptive_table = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
        loc, scale = np.mean(adaptive_table), np.std(adaptive_table)
        x = filtered_table_omega(dico_omega, nearly_neutral_dico, gene)[:, 1]
        assert (not np.isnan(x).any())
        weigths_nearly_neutral = norm.pdf(x, loc=loc, scale=scale)
        weigths_nearly_neutral /= norm.pdf(x, loc=np.mean(x), scale=np.std(x))
        weigths_nearly_neutral /= np.sum(weigths_nearly_neutral)

    for rep in range(1, args.rep + 1):
        if gene:
            adaptive_set = subsample_genes(adaptive_dico, args.nbr_genes, None)
            nearly_neutral_set = subsample_genes(nearly_neutral_dico, args.nbr_genes, weigths_nearly_neutral)
        else:
            adaptive_set = subsample_sites(adaptive_dico, args.nbr_sites, None)
            nearly_neutral_set = subsample_sites(nearly_neutral_dico, args.nbr_sites, weigths_nearly_neutral)

        w = filtered_table_omega(dico_omega, nearly_neutral_set, gene)[:, 1]
        w_0 = filtered_table_omega(dico_omega_0, nearly_neutral_set, gene)[:, 1]
        w_A = w - w_0
        f.write(txt.format("", np.mean(w), np.std(w), "resampled nearly-neutral", args.granularity))
        f.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled nearly-neutral", args.granularity))
        f.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled nearly-neutral", args.granularity))

        w = filtered_table_omega(dico_omega, adaptive_set, gene)[:, 1]
        w_0 = filtered_table_omega(dico_omega_0, adaptive_set, gene)[:, 1]
        w_A = w - w_0
        f.write(txt.format("", np.mean(w), np.std(w), "resampled adaptive", args.granularity))
        f.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled adaptive", args.granularity))
        f.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled adaptive", args.granularity))

        dfe_alpha("{0}_ADAPTIVE_{1}.txt".format(args.output, rep), df_snp, sample_size, adaptive_set, gene,
                  args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder, args.dfem_path)

        dfe_alpha("{0}_NEARLY_NEUTRAL_{1}.txt".format(args.output, rep), df_snp, sample_size, nearly_neutral_set, gene,
                  args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder, args.dfem_path)

    f.close()
