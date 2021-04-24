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
    parser.add_argument('--nbr_bins', required=True, type=int, dest="nbr_bins", help="Number of bins")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")

    args = parser.parse_args()

    gene = args.granularity.lower() == "gene"
    dico_omega_0, dico_omega = build_divergence_dico(args.div_folder, gene_level=gene)
    dico_alignments = load_alignments(dico_omega, args.focal_species, args.sister_species, args.ali_folder)
    dico_bins = bin_dataset(dico_omega_0, dico_omega, args.nbr_bins, gene_level=gene, filter_set=set(dico_alignments))

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf)
    np.random.seed(args.seed)

    txt = "E[w{0}]={1:.3g}, Var[w{0}]={2:.3g} for bin {3} ({4} level) set\n"
    f = open(args.output + "bins.txt", 'w')
    for bin_id, dico in enumerate(dico_bins):
        w = filtered_table_omega(dico_omega, dico, gene)[:, 1]
        w_0 = filtered_table_omega(dico_omega_0, dico, gene)[:, 1]
        w_A = w - w_0
        f.write(txt.format("", np.mean(w), np.std(w), bin_id, args.granularity))
        f.write(txt.format("0", np.mean(w_0), np.std(w_0), bin_id, args.granularity))
        f.write(txt.format("A", np.mean(w_A), np.std(w_A), bin_id, args.granularity))

    for rep in range(1, args.rep + 1):
        for bin_id, dico in enumerate(dico_bins):
            if gene:
                resampled_set = subsample_genes(dico, args.nbr_genes, None)
            else:
                resampled_set = subsample_sites(dico, args.nbr_sites, None)

            w = filtered_table_omega(dico_omega, resampled_set, gene)[:, 1]
            w_0 = filtered_table_omega(dico_omega_0, resampled_set, gene)[:, 1]
            w_A = w - w_0
            f.write(txt.format("", np.mean(w), np.std(w), "resampled ", args.granularity))
            f.write(txt.format("0", np.mean(w_0), np.std(w_0), "resampled nearly-neutral", args.granularity))
            f.write(txt.format("A", np.mean(w_A), np.std(w_A), "resampled nearly-neutral", args.granularity))

            dfe_alpha("{0}{1}_{2}.txt".format(args.output, bin_id + 1, rep), df_snp, sample_size, resampled_set, gene,
                      args.focal_species, args.sister_species, dico_alignments, fix_poly, args.tmp_folder,
                      args.dfem_path)
    f.close()