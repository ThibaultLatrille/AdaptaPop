import argparse
import bz2
import _pickle as cPickle
from libraries import *
from scipy.stats import norm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--level', required=False, type=str, default="gene", dest="level",
                        help="Gene or site level")
    parser.add_argument('-m', '--method', required=False, type=str, default='MutSel', dest="method",
                        help="Method to detect adaptation (MutSel, Classical, MutSelExclu)")
    parser.add_argument('--nbr_sites', required=False, type=int, default=-1, dest="nbr_sites",
                        help="Number of sites to subsample")
    parser.add_argument('--nbr_genes', required=False, type=int, default=-1, dest="nbr_genes",
                        help="Number of genes to subsample")
    parser.add_argument('--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('--weighted', required=False, type=str, default='False', dest="weighted",
                        help="Weighted sampling such as to control for omega")
    parser.add_argument('--seed', required=False, type=int, default=123456789, dest="seed",
                        help="Seed for random generator")
    parser.add_argument('--pickle', required=True, type=str, dest="pickle", help="Pickle file")
    parser.add_argument('--dfe_path', required=False, type=str, dest="dfe_path", nargs="+", help="Executable path")
    parser.add_argument('--sfs', required=False, type=str, dest="sfs", default="folded", help="unfolded or folded")
    parser.add_argument('--subsample', required=False, type=int, default=-1, dest="subsample", help="Subsample SFS")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    gene = "gene" in args.level.lower()
    assert args.method in ["MutSel", "Classical", "MutSelExclu"]

    polarize_snps = args.sfs == "unfolded"
    pickle_file = bz2.BZ2File(args.pickle, 'rb')
    filter_set, dico_omega_0, dico_omega, dico_alignments = cPickle.load(pickle_file)
    pickle_file.close()
    print("Data loaded")
    adaptive_dico, nn_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=gene, filter_set=filter_set,
                                               method=args.method)
    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf, polarize_snps)
    adaptive_dico = {k: v for k, v in adaptive_dico.items() if k in df_snp.groups}
    nn_dico = {k: v for k, v in nn_dico.items() if k in df_snp.groups}

    print("Data classified")
    if gene:
        if args.nbr_genes == -1 or args.nbr_genes > len(adaptive_dico):
            args.nbr_genes = len(adaptive_dico)
        print('{0} adaptive genes'.format(len(adaptive_dico)))
        print('{0} nearly-neutral genes'.format(len(nn_dico)))
    else:
        adapta_sites = sum([len(v) for v in adaptive_dico.values()])
        if args.nbr_sites == -1 or args.nbr_sites > adapta_sites:
            args.nbr_sites = adapta_sites
        print('{0} adaptive sites'.format(adapta_sites))
        print('{0} nearly-neutral sites'.format(sum([len(v) for v in nn_dico.values()])))

    np.random.seed(args.seed)

    weigths_nearly_neutral = None
    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for {3} {4} set\n"
    f_sample = open(args.output + ".txt", 'w')
    w = filtered_table_omega(dico_omega, nn_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, nn_dico, gene)[:, 1]
    w_A = w - w_0
    f_sample.write(txt.format("", np.mean(w), np.std(w), "all nearly-neutral", args.level))
    f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "all nearly-neutral", args.level))
    f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "all nearly-neutral", args.level))

    w = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, adaptive_dico, gene)[:, 1]
    w_A = w - w_0
    f_sample.write(txt.format("", np.mean(w), np.std(w), "all adaptive", args.level))
    f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "all adaptive", args.level))
    f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "all adaptive", args.level))

    errors = gzip.open(args.tmp_folder + "/errors.txt.gz", 'wt')
    if args.output.split("/")[-1] == "1":
        print("{0}_ADAPTIVE".format(args.output))
        dfe_alpha("{0}_ADAPTIVE".format(args.output), df_snp, sample_size, args.subsample, adaptive_dico, gene,
                  args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder,
                  args.dfe_path, polarize_snps, errors)

    if args.weighted.lower() == "true":
        assert (args.method == "MutSelExclu")
        print("Weighted sampling controlling for w")
        adaptive_table = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
        loc, scale = np.mean(adaptive_table), np.std(adaptive_table)
        x = filtered_table_omega(dico_omega, nn_dico, gene)[:, 1]
        assert (not np.isnan(x).any())
        weigths_nearly_neutral = norm.pdf(x, loc=loc, scale=scale)
        weigths_nearly_neutral /= norm.pdf(x, loc=np.mean(x), scale=np.std(x))
        weigths_nearly_neutral /= np.sum(weigths_nearly_neutral)

    prefix = "{0}_NEARLY_NEUTRAL".format(args.output)
    rep, nbr_errors = 1, 1
    while rep < args.rep + 1 and nbr_errors < args.rep + 1:
        print("{0}_{1}".format(prefix, rep))
        if gene:
            subset = subsample_genes(nn_dico, args.nbr_genes, weigths_nearly_neutral, replace=False)
        else:
            subset = subsample_sites(nn_dico, args.nbr_sites, weigths_nearly_neutral, replace=False)

        w = filtered_table_omega(dico_omega, subset, gene)[:, 1]
        w_0 = filtered_table_omega(dico_omega_0, subset, gene)[:, 1]
        w_A = w - w_0
        f_sample.write(txt.format("", np.mean(w), np.std(w), "resampled nearly-neutral", args.level))
        f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled nearly-neutral", args.level))
        f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled nearly-neutral", args.level))

        if dfe_alpha("{0}_{1}".format(prefix, rep), df_snp, sample_size, args.subsample, subset, gene,
                     args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder,
                     args.dfe_path, polarize_snps, errors):
            rep += 1
        else:
            nbr_errors += 1

    errors.close()
    f_sample.close()
