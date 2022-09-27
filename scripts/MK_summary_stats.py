import argparse
import bz2
import _pickle as cPickle

from libraries import *


def MK_stats(df, k_tot, k, ensg, sp_focal, sp_sister, ali_dico, fixed_poly, folder, filepath, error_f):
    p_non_syn, p_syn = 0, 0
    s1, s2 = [], []
    ensg_fixedpoly = dict(fixed_poly[ensg])

    dff = df.get_group(ensg)
    for row in dff.itertuples(index=False):
        if row.TYPE == "Stop":
            continue

        if k_tot > k:
            daf = np.random.hypergeometric(row.COUNT, k_tot - row.COUNT, k)
            if (daf == 0) or (daf == k):
                continue

        if row.TYPE == "Syn":
            p_syn += 1
        else:
            assert row.TYPE == "NonSyn"
            p_non_syn += 1

    seq_dico = filter_positions(ali_dico[ensg], sp_focal, ensg_fixedpoly, True, None)
    seq1 = SeqIO.SeqRecord(id=sp_focal, name="", description="", seq=Seq.Seq(seq_dico[sp_focal]))
    seq2 = SeqIO.SeqRecord(id=sp_sister, name="", description="", seq=Seq.Seq(seq_dico[sp_sister]))

    if seq1.seq == seq2.seq:
        error_f.write("ER1: identical sequences. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    l_non_syn, d_non_syn, l_syn, d_syn = run_yn00(seq1, seq2, folder, filepath)
    if l_non_syn == d_non_syn == l_syn == d_syn == 0:
        print('error in Yn00')
        error_f.write("ER2: Yn00 failed. Sequences are:\n" + "".join(s1) + "\n" + "".join(s2) + "\n")
        return False

    return l_non_syn, d_non_syn, p_non_syn, l_syn, d_syn, p_syn


def main(args):
    os.makedirs(args.tmp_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    pickle_file = bz2.BZ2File(args.pickle, 'rb')
    filter_set, dico_omega_0, dico_omega, dico_alignments, dico_trid = cPickle.load(pickle_file)
    annot_dico = {ensg: [trid, chrom, strand] for ensg, trid, chrom, strand in
                  zip(dico_trid['ENSG'], dico_trid['TRID'], dico_trid['CHR'], dico_trid['STRAND'])}
    pickle_file.close()

    df_snp, fix_poly, sample_size = snp_data_frame(args.vcf, False)
    print(f"Number of genes in the filter set: {len(filter_set)}")
    filter_set = list(sorted(set(df_snp.groups) & filter_set))
    print(f"Number of genes in the VCF: {len(df_snp.groups)}")
    print(f"Number of genes in the intersection: {len(filter_set)}")

    for ensg in filter_set:
        assert ensg.count("_") == 1
    df_out = pd.DataFrame(columns=["ENSG", "NAME", "TRID", "CHR", "STRAND", "L_non_syn", "D_non_syn", "P_non_syn",
                                   "L_syn", "D_syn", "P_syn"])
    # ENSG is the gene ID on Ensembl shared by all species (in the file name of OrthoMam alignment)
    # NAME is the gene name  shared by all species (in the file name of OrthoMam alignment)
    # TRID is the transcript ID of the gene, specific to the focal species (found in the .xml files of OrthoMam)
    # CHR is the chromosome on which the gene is located
    # STRAND is the strand on which the gene is located (+ if the same as the reference genome, - otherwise)
    # L_non_syn is the number of non-synonymous sites on which the substitutions and polymorphisms are called
    # D_non_syn is the number of non-synonymous substitutions (can be 0)
    # P_non_syn is the number of non-synonymous polymorphisms (can be 0)
    # L_syn is the number of synonymous sites on which the substitutions and polymorphisms are called
    # D_syn is the number of synonymous substitutions (can be 0)
    # P_syn is the number of synonymous polymorphisms (can be 0)
    # dN is D_non_syn / L_non_syn, dS is D_syn / L_syn, and dN/dS is thus D_non_syn / D_syn
    # pN is P_non_syn / L_non_syn, pS is P_syn / L_syn, and pN/pS is thus P_non_syn / P_syn
    pct = 0
    for ensg_index, ensg in enumerate(filter_set):
        new_pct = int(100 * (ensg_index + 1) / len(filter_set))
        if pct != new_pct:
            pct = new_pct
            print(f"{new_pct}%")
        error_f = open(f'{args.tmp_folder}errors.txt', 'w')
        filepath = f'{args.tmp_folder}{ensg}'
        stat = MK_stats(
            df_snp, sample_size, args.subsample, ensg, args.focal_species, args.sister_species, dico_alignments,
            fix_poly, args.tmp_folder, filepath, error_f)
        if stat is False:
            continue
        # l_non_syn, d_non_syn, p_non_syn, l_syn, d_syn, p_syn = stat
        # dn, ds = d_non_syn / l_non_syn, d_syn / l_syn
        # dnds = dn / ds if ds != 0 else np.nan
        # pn, ps = p_non_syn / l_non_syn, p_syn / l_syn
        # pnps = pn / ps if ps != 0 else np.nan
        trid, chrom, strand = annot_dico[ensg]
        df_out.loc[len(df_out)] = ensg.split("_") + [trid, chrom, strand] + list(stat)  # + [dnds, pnps]
    df_out.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', required=True, type=str, dest="vcf", help="VCF annotated file")
    parser.add_argument('--focal_species', required=True, type=str, dest="focal_species",
                        help="Focal species (containing polymorphism)")
    parser.add_argument('--sister_species', required=True, type=str, dest="sister_species",
                        help="Sister species (not containing polymorphism)")
    parser.add_argument('--pickle', required=True, type=str, dest="pickle", help="Pickle file")
    parser.add_argument('--subsample', required=False, type=int, default=-1, dest="subsample", help="Subsample SFS")
    parser.add_argument('--tmp_folder', required=True, type=str, dest="tmp_folder", help="Temporary files path")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output path")
    main(parser.parse_args())
