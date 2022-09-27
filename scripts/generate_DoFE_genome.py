import argparse
import bz2
import _pickle as cPickle
from libraries import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('--seed', required=False, type=int, default=123456789, dest="seed",
                        help="Seed for random generator")
    parser.add_argument('--pickle', required=True, type=str, dest="pickle", help="Pickle file")
    parser.add_argument('--dfe_path', required=False, type=str, dest="dfe_path", nargs="+", help="Executable path")
    parser.add_argument('--sfs', required=False, type=str, dest="sfs", default="folded", help="unfolded or folded")
    parser.add_argument('--subsample', required=False, type=int, default=-1, dest="subsample", help="Subsample SFS")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    polarize_snps = args.sfs == "unfolded"
    pickle_file = bz2.BZ2File(args.pickle, 'rb')
    filter_set, dico_omega_0, dico_omega, dico_alignments, dico_trid = cPickle.load(pickle_file)
    pickle_file.close()
    print("Data loaded")

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf, polarize_snps)
    filter_set = set(df_snp.groups) & filter_set
    np.random.seed(args.seed)

    gene_len_array = np.array([len(dico_alignments[ensg][args.focal_species]) for ensg in filter_set])
    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for all genome\n"
    f_sample = open(args.output + ".txt", 'w')

    w_array = filtered_table_omega(dico_omega, filter_set, True)[:, 1].T
    assert w_array.shape == gene_len_array.shape
    w = (gene_len_array * w_array) / np.sum(gene_len_array)

    w_0_array = filtered_table_omega(dico_omega_0, filter_set, True)[:, 1].T
    assert w_0_array.shape == gene_len_array.shape
    w_0 = (gene_len_array * w_0_array) / np.sum(gene_len_array)
    w_A = w - w_0
    f_sample.write(txt.format("", np.mean(w), np.std(w)))
    f_sample.write(txt.format("0", np.mean(w_0), np.std(w_0)))
    f_sample.write(txt.format("A", np.mean(w_A), np.std(w_A)))

    errors = gzip.open(args.tmp_folder + "/errors.txt.gz", 'wt')
    dico_ensg = {i: None for i in filter_set}
    dfe_alpha(args.output, df_snp, sample_size, args.subsample, dico_ensg, True,
              args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder,
              args.dfe_path, polarize_snps, errors)

    errors.close()
    f_sample.close()
