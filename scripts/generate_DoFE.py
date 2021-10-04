import argparse
import bz2
import _pickle as cPickle
from libraries import *
from scipy.stats import norm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--granularity', required=False, type=str, default="gene", dest="granularity",
                        help="Gene or site level")
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
    parser.add_argument('--dfe_path', required=True, type=str, dest="dfe_path", nargs="+", help="Executable path")
    parser.add_argument('--sfs', required=False, type=str, dest="sfs", default="folded", help="unfolded or folded")
    parser.add_argument('--gather', required=False, type=str, dest="gather", default="false", help="Gather SFS")
    parser.add_argument('--subsample', required=False, type=int, default=-1, dest="subsample", help="Subsample SFS")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    gather_sfs = args.gather.lower() == "true"
    is_unfolded = args.sfs == "unfolded"
    pickle_file = bz2.BZ2File(args.pickle, 'rb')
    filter_set, dico_omega_0, dico_omega, dico_alignments = cPickle.load(pickle_file)
    pickle_file.close()
    print("Data loaded")
    _, adaptive_dico, epistasis_dico, nearly_neutral_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=gene,
                                                                              filter_set=filter_set)
    print("Data classified")
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

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf, is_unfolded, args.subsample * 2)
    np.random.seed(args.seed)

    weigths_nearly_neutral = None
    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for {3} {4} set\n"
    f_sample = open(args.output + ".txt", 'w')
    w = filtered_table_omega(dico_omega, nearly_neutral_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, nearly_neutral_dico, gene)[:, 1]
    w_A = w - w_0
    f_sample.write(txt.format("", np.mean(w), np.std(w), "all nearly-neutral", args.granularity))
    f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "all nearly-neutral", args.granularity))
    f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "all nearly-neutral", args.granularity))

    w = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
    w_0 = filtered_table_omega(dico_omega_0, adaptive_dico, gene)[:, 1]
    w_A = w - w_0
    f_sample.write(txt.format("", np.mean(w), np.std(w), "all adaptive", args.granularity))
    f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "all adaptive", args.granularity))
    f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "all adaptive", args.granularity))

    errors = gzip.open(args.tmp_folder + "/errors.txt.gz", 'wt')
    if args.weighted.lower() == "true":
        print("Weighted sampling controlling for w")
        adaptive_table = filtered_table_omega(dico_omega, adaptive_dico, gene)[:, 1]
        loc, scale = np.mean(adaptive_table), np.std(adaptive_table)
        x = filtered_table_omega(dico_omega, nearly_neutral_dico, gene)[:, 1]
        assert (not np.isnan(x).any())
        weigths_nearly_neutral = norm.pdf(x, loc=loc, scale=scale)
        weigths_nearly_neutral /= norm.pdf(x, loc=np.mean(x), scale=np.std(x))
        weigths_nearly_neutral /= np.sum(weigths_nearly_neutral)

    for adapt in [True, False]:
        txt_set = "ADAPTIVE" if adapt else "NEARLY_NEUTRAL"
        prefix = "{0}_{1}".format(args.output, txt_set)
        sfs_polyDFE = {"SFSn": "", "SFSs": ""}
        rep = 1

        while rep < args.rep + 1:
            print("{0}_{1}".format(prefix, rep))
            if gene:
                if adapt:
                    subset = subsample_genes(adaptive_dico, args.nbr_genes, None, replace=True)
                else:
                    subset = subsample_genes(nearly_neutral_dico, args.nbr_genes, weigths_nearly_neutral, replace=False)
            else:
                if adapt:
                    subset = subsample_sites(adaptive_dico, args.nbr_sites, None, replace=True)
                else:
                    subset = subsample_sites(nearly_neutral_dico, args.nbr_sites, weigths_nearly_neutral, replace=False)

            w = filtered_table_omega(dico_omega, subset, gene)[:, 1]
            w_0 = filtered_table_omega(dico_omega_0, subset, gene)[:, 1]
            w_A = w - w_0
            f_sample.write(txt.format("", np.mean(w), np.std(w), "resampled " + txt_set.lower(), args.granularity))
            f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled " + txt_set.lower(), args.granularity))
            f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled " + txt_set.lower(), args.granularity))

            if dfe_alpha("{0}_{1}".format(prefix, rep), df_snp, sample_size, subset, gene,
                         args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder,
                         args.dfe_path, is_unfolded, sfs_polyDFE, errors, gather_sfs):
                rep += 1

        if sfs_polyDFE["SFSs"] != "" and gather_sfs:
            sfs_file = open(prefix + ".sfs", 'w')
            sfs_file.write("#{0}+{1}".format(args.focal_species, args.sister_species) + "\n")
            sfs_file.write("{0} {0} {1}".format(args.rep, sample_size) + "\n")
            sfs_file.write(sfs_polyDFE["SFSs"])
            sfs_file.write(sfs_polyDFE["SFSn"])
            sfs_file.close()
            dfe_path = set([p for p in args.dfe_path if "polyDFE" in p]).pop()
            os.system(polyDFE_cmd.format(dfe_path, prefix, prefix + "_polyDFE"))

    errors.close()
    f_sample.close()