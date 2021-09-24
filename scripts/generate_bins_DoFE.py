import argparse
from libraries import *
import bz2
import _pickle as cPickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--granularity', required=False, type=str, default="gene", dest="granularity",
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
    parser.add_argument('--pickle', required=True, type=str, dest="pickle", help="Annotation of ENSG")
    parser.add_argument('--dfe_path', required=True, type=str, dest="dfe_path", nargs="+", help="Executable path")
    parser.add_argument('--subsample', required=False, type=int, default=-1, dest="subsample",  help="Subsample SFS")
    parser.add_argument('--sfs', required=False, type=str, dest="sfs", default="folded", help="unfolded or folded")
    parser.add_argument('--gather', required=False, type=str, dest="gather", default="false", help="Gather SFS")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--rep', required=True, type=int, dest="rep", help="Number of replicates")
    parser.add_argument('--nbr_bins', required=True, type=int, dest="nbr_bins", help="Number of bins")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    gather_sfs = args.gather.lower() == "true"
    is_unfolded = args.sfs == "unfolded"
    if args.nbr_sites <= 0: args.nbr_sites = float('inf')
    if args.nbr_genes <= 0: args.nbr_genes = float('inf')
    pickle_file = bz2.BZ2File(args.pickle, 'rb')
    filter_set, dico_omega_0, dico_omega, dico_alignments = cPickle.load(pickle_file)
    pickle_file.close()
    print("Data loaded")
    _, _, _, nn_dico, _ = split_outliers(dico_omega_0, dico_omega, gene_level=gene, filter_set=filter_set)
    dico_bins = bin_dataset(nn_dico, dico_omega_0, args.nbr_bins, gene_level=gene)
    print("Data binned")

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf, is_unfolded, args.subsample * 2)
    np.random.seed(args.seed)

    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for bin {3} ({4} level) set\n"
    f_sample = open(args.output + "bins.txt", 'w')
    errors = gzip.open(args.tmp_folder + "/errors.txt.gz", 'wt')

    for bin_id, dico in enumerate(dico_bins):
        w = filtered_table_omega(dico_omega, dico, gene)[:, 1]
        w_0 = filtered_table_omega(dico_omega_0, dico, gene)[:, 1]
        w_A = w - w_0
        f_sample.write(txt.format("", np.mean(w), np.std(w), bin_id, args.granularity))
        f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), bin_id, args.granularity))
        f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), bin_id, args.granularity))

        prefix = "{0}{1}".format(args.output, np.mean(w_0))
        sfs_polyDFE = {"SFSn": "", "SFSs": ""}
        rep = 1
        while rep < args.rep + 1:
            print("{0}_{1}".format(prefix, rep))
            if gene:
                resampled_set = subsample_genes(dico, min(len(w_0), args.nbr_genes), None, replace=True)
            else:
                resampled_set = subsample_sites(dico, min(len(w_0), args.nbr_sites), None, replace=True)

            w = filtered_table_omega(dico_omega, resampled_set, gene)[:, 1]
            w_0 = filtered_table_omega(dico_omega_0, resampled_set, gene)[:, 1]
            w_A = w - w_0
            f_sample.write(txt.format("", np.mean(w), np.std(w), "resampled ", args.granularity))
            f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled nearly-neutral", args.granularity))
            f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled nearly-neutral", args.granularity))

            if dfe_alpha("{0}_{1}".format(prefix, rep), df_snp, sample_size, resampled_set, gene,
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
