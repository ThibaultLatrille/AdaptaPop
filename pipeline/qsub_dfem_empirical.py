import os
import argparse
import numpy as np
from math import floor
from collections import Counter
from cds_libraries import load_pb_table, load_snp_table
import scipy.stats as st


def sample_sfs(sfs, n_sample, n, replacement, theoretical=True):
    assert n_sample <= n
    fold = int(floor(n_sample // 2))
    if theoretical:
        out = {i: 0. for i in range(n_sample+1)}
        for snp in sfs:
            x = np.arange(0, min(snp, n_sample) + 1)
            pdf = st.hypergeom(n, snp, n_sample).pmf(x)
            for k in x:
                if k < fold + 1:
                    out[k] += pdf[k]
                else:
                    out[n_sample-k] += pdf[k]
        assert sum(list(out.values())[fold+1:]) == 0
        return out
    else:
        if replacement:
            sampled_sfs = np.random.binomial(n_sample, [i/n for i in sfs])
        else:
            sampled_sfs = np.random.hypergeometric(sfs, [n - i for i in sfs], n_sample)
        return Counter([i if i < fold + 1 else n_sample - i for i in list(sampled_sfs) if i > 0])


def dfem(group_mk, dfem_path, sample_nbr, sfs_filestr):
    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = [], []
    nbr_alleles = 5008
    _sample_size = nbr_alleles if nbr_alleles < sample_nbr else sample_nbr
    if _sample_size % 2 == 1:
        _sample_size -= 1
    replacement = False
    for mk in group_mk:
        sites_s += mk["SitesS"]
        sites_n += mk["SitesN"]
        ds += mk["Ds"]
        dn += mk["Dn"]
        sfs_s.extend(mk["SFSs"])
        sfs_n.extend(mk["SFSn"])
        ps += mk["Ps"]
        pn += mk["Pn"]

    sfs_s_sample = sample_sfs(sfs_s, _sample_size, nbr_alleles, replacement)
    sfs_n_sample = sample_sfs(sfs_n, _sample_size, nbr_alleles, replacement)

    sfs_list = ["Summed", str(_sample_size)]
    for sfs, nbr_site in [(sfs_n_sample, sites_n), (sfs_s_sample, sites_s)]:
        range_sfs = range(1, int(floor(_sample_size // 2)) + 1)
        assert len(range_sfs) * 2 == _sample_size
        sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]

    sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
    sfs_file = open("{0}/{1}".format(dfem_path, sfs_filestr), 'w')
    sfs_file.write("Homo\n" + "\t".join(sfs_list) + "\n")
    sfs_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--name', required=True, type=int,
                        dest="n", metavar="<qsub id>",
                        help="Qsub id")
    args = parser.parse_args()
    qsub_name = "empirical_{0}".format(args.n)
    qsub_id = args.n

    data_path = "/panhome/tlatrill/AdaptaPop/data"
    version = "88"
    GRCh = "GRCh38"

    snp_table = load_snp_table(data_path, version, GRCh)
    pb_table = load_pb_table(data_path)
    pb_cds_ids = set(cds_id[2:-1] for cds_id in pb_table["CdsId"])

    outlier_set = set()
    cds_file = open('{0}/outliers.out'.format(data_path), 'r')
    cds_file.readline()
    cds_file.readline()
    for line in cds_file:
        if line != '\n':
            outlier_set.add(line.split("\t")[0])
    cds_file.close()

    snp_set = set([k for k, v in snp_table.items() if v["Chromosome"] not in ["X", "Y", "MT"]])
    outlier_set = outlier_set & snp_set
    non_outliers_list = list((pb_cds_ids - outlier_set) & snp_set)
    print('{0} outliers'.format(len(outlier_set)))
    print('{0} non outliers'.format(len(non_outliers_list)))

    model = "GammaExpo"
    dfe_path = "/panhome/tlatrill/AdaptaPop/dfem"

    sample_size = 12
    nbr_repeats = 500

    dfem_non_outliers = []
    for repeat_id in range(nbr_repeats):
        group = [snp_table[cds_id] for cds_id in np.random.choice(non_outliers_list, len(outlier_set), replace=False)]

        sfs_filename = "sfs_{0}_{1}.txt".format(qsub_name, repeat_id)
        out_filename = "out_{0}_{1}.csv".format(qsub_name, repeat_id)
        dfem(group, dfe_path, sample_size, sfs_filename)

        os.system("{0}/dfem  -in {0}/{1} -out {0}/{2} -model {3}\n".format(dfe_path, sfs_filename, out_filename, model))

        out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
        headline = out_file.readline().split(',')
        index_model = headline.index("model")

        for line in out_file:
            line_split = line.split(',')
            if line_split[index_model] == model:
                omega_a = float(line_split[headline.index("omegaA")])
                dfem_non_outliers.append(omega_a)
        out_file.close()
        os.remove("{0}/{1}".format(dfe_path, out_filename))
        os.remove("{0}/{1}".format(dfe_path, sfs_filename))

    file = open("{0}/dfem_non_outliers_{0}.p".format(data_path, qsub_name), 'wb')
    pickle.dump(dfem_non_outliers, file)
    file.close()
    print(qsub_name)
