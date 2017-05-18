import numpy as np
from math import floor
from collections import Counter
from cds_libraries import load_pb_table, load_snp_table
import scipy.stats as st

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


def dfem(group, dfe_path, sample_size, sfs_filename):
    sites_n, sites_s, dn, ds, pn, ps = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = [], []
    nbr_alleles = 5008
    sample_size = nbr_alleles if nbr_alleles < sample_size else sample_size
    if sample_size % 2 == 1:
        sample_size -= 1
    replacement = False
    for mk in group:
        sites_s += mk["SitesS"]
        sites_n += mk["SitesN"]
        ds += mk["Ds"]
        dn += mk["Dn"]
        sfs_s.extend(mk["SFSs"])
        sfs_n.extend(mk["SFSn"])
        ps += mk["Ps"]
        pn += mk["Pn"]

    sfs_s_sample = sample_sfs(sfs_s, sample_size, nbr_alleles, replacement)
    sfs_n_sample = sample_sfs(sfs_n, sample_size, nbr_alleles, replacement)

    sfs_list = ["Summed", str(sample_size)]
    for sfs, nbr_site in [(sfs_n_sample, sites_n), (sfs_s_sample, sites_s)]:
        range_sfs = range(1, int(floor(sample_size // 2)) + 1)
        assert len(range_sfs) * 2 == sample_size
        sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]

    sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
    sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
    sfs_file.write("Homo\n")
    sfs_file.write("\t".join(sfs_list) + "\n")
    sfs_file.close()


model = "GammaExpo"
dfe_path = "/panhome/tlatrill/AdaptaPop/dfem"

sample_size = 12
nbr_qsubs = 200
nbr_repeats = 2000

for qsub_id in range(nbr_qsubs):
    qsub_name = "empirical_{0}".format(qsub_id)
    qsub = open("{0}/qsub/{1}.pbs".format(dfe_path, qsub_name), 'w')
    qsub.writelines("#!/bin/bash\n")
    qsub.writelines("#\n")
    qsub.writelines("#PBS -q q1day\n")
    qsub.writelines("#PBS -l nodes=1:ppn=1,mem=1gb\n")
    qsub.writelines("#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(qsub_name))
    qsub.writelines("#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(qsub_name))
    qsub.writelines("#PBS -j oe\n")
    qsub.writelines("#PBS -W umask=022\n")
    qsub.writelines("#PBS -r n\n")
    qsub.writelines("#PBS -r n\n")
    for repeat_id in range(nbr_repeats):
        group = [snp_table[cds_id] for cds_id in np.random.choice(non_outliers_list, len(outlier_set), replace=False)]

        sfs_filename = "sfs_sample_{0}_{1}.txt".format(qsub_name, repeat_id)
        out_filename = "out_sample_{0}_{1}.csv".format(qsub_name, repeat_id)
        dfem(group, dfe_path, sample_size, sfs_filename)

        qsub.writelines("{0}/dfem  -in {0}/{1} -out {0}/{2} -model {3}\n".format(dfe_path, sfs_filename, out_filename, model))
    qsub.writelines("rm {0}/qsub/{1}.pbs\n".format(dfe_path, qsub_name))
    qsub.close()
    print(qsub_name)


nbr_repeats = 200
groups = [[snp_table[cds_id] for cds_id in np.random.choice(non_outliers_list, len(outlier_set), replace=False)] for _ in range(nbr_repeats)]
sample_size_list = [int(i) for i in np.logspace(np.log10(12), np.log10(192), 30)]

for qsub_id, group in enumerate(groups):
    qsub_name = "subsampling_{0}".format(qsub_id)
    qsub = open("{0}/qsub/{1}.pbs".format(dfe_path, qsub_name), 'w')
    qsub.writelines("#!/bin/bash\n")
    qsub.writelines("#\n")
    qsub.writelines("#PBS -q q1day\n")
    qsub.writelines("#PBS -l nodes=1:ppn=1,mem=1gb\n")
    qsub.writelines("#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(qsub_name))
    qsub.writelines("#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(qsub_name))
    qsub.writelines("#PBS -j oe\n")
    qsub.writelines("#PBS -W umask=022\n")
    qsub.writelines("#PBS -r n\n")
    qsub.writelines("#PBS -r n\n")
    for sample_size in sample_size_list:
        sfs_filename = "sfs_{0}_{1}.txt".format(qsub_name, sample_size)
        out_filename = "out_{0}_{1}.csv".format(qsub_name, sample_size)
        dfem(group, dfe_path, sample_size, sfs_filename)

        qsub.writelines("{0}/dfem  -in {0}/{1} -out {0}/{2} -model {3}\n".format(dfe_path, sfs_filename, out_filename, model))
    qsub.writelines("rm {0}/qsub/{1}.pbs\n".format(dfe_path, qsub_name))
    qsub.close()
    print(qsub_name)
