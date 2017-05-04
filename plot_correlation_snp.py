import os
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
from math import floor
import numpy as np
from collections import Counter
from cds_libraries import load_snp_table, load_pb_table

data_path = "/mnt/sda1/AdaptaPop/data"
phase3 = True
version = "88"
GRCh = "GRCh38"
id_str = "jk_intra_1"

pb_table = load_pb_table(data_path)
pb_cds_ids = set(cds_id[2:-1] for cds_id in pb_table["CdsId"])
keep_cds_probability = 1


def sample_sfs(sfs, n_sample, n, replacement):
    assert n_sample <= n
    if replacement:
        sampled_sfs = np.random.binomial(n_sample, [i/n for i in sfs])
    else:
        sampled_sfs = np.random.hypergeometric(sfs, [n - i for i in sfs], n_sample)
    return [i if i < n_sample // 2 + 1 else n_sample // 2 - i for i in list(sampled_sfs) if i > 0]


def dfe_alpha(group, sample_size):
    sites_n, sites_s, dn, ds, pn, ps, sfs_n_len, sfs_s_len, = 0, 0, 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = Counter(), Counter()
    nbr_alleles = 5008 if phase3 else 2178
    sample_size = nbr_alleles if nbr_alleles < sample_size else sample_size
    if sample_size % 2 == 1:
        sample_size -= 1
    replacement = False
    for mk in group:
        if mk["Chromosome"] not in ["X", "Y", "MT"]:
            sites_s += mk["SitesS"]
            sites_n += mk["SitesN"]
            ds += mk["Ds"]
            dn += mk["Dn"]
            sfs_s_sample = sample_sfs(mk["SFSs"], sample_size, nbr_alleles, replacement)
            sfs_n_sample = sample_sfs(mk["SFSn"], sample_size, nbr_alleles, replacement)
            if sample_size == nbr_alleles and not replacement:
                assert len(sfs_s_sample) == len(mk["SFSs"])
                assert sum(sfs_s_sample) == sum(mk["SFSs"])
                for i, j in zip(sfs_s_sample, mk["SFSs"]):
                    assert i == j
            sfs_s.update(sfs_s_sample)
            sfs_n.update(sfs_n_sample)
            sfs_s_len += len(sfs_s_sample)
            sfs_n_len += len(sfs_n_sample)
            ps += mk["Ps"]
            pn += mk["Pn"]
    dfe_path = "/home/thibault/Tools/dfem"
    sfs_list = ["Summed", str(sample_size)]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        range_sfs = range(1, int(floor(sample_size // 2))+1)
        assert len(range_sfs) * 2 == sample_size
        sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]
    sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
    sfs_filename = "sfs_{0}.txt".format(id_str)
    sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
    sfs_file.write("Homo\n")
    sfs_file.write("\t".join(sfs_list) + "\n")
    sfs_file.close()

    out_filename = "output_{0}.csv".format(id_str)
    os.system("{0}/dfem  -in {0}/{1} -out {0}/{2} -model GammaExpo".format(dfe_path, sfs_filename, out_filename))
    out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
    headline = out_file.readline().split(',')
    index_alpha = headline.index("omegaA")
    index_omega = headline.index("alpha")
    index_model = headline.index("model")
    alpha_list = []
    omega_list = []
    for line in out_file:
        line_split = line.split(',')
        if line_split[index_model] == "GammaExpo":
            alpha_str = line_split[index_alpha]
            if alpha_str != "NA":
                alpha_list.append(float(alpha_str))
            omega_str = line_split[index_omega]
            if omega_str != "NA":
                omega_list.append(float(omega_str))
    out_file.close()
    return sum(alpha_list) / len(alpha_list), sum(omega_list) / len(omega_list), \
           pn * sites_s / (ps * sites_n), dn * sites_s / (ds * sites_n)

data_path = "/mnt/sda1/AdaptaPop/data"
my_dpi = 96
plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
snp_table = load_snp_table(data_path, version, GRCh, 1)


sample_size_list = [64]*8
cds_list = [mk for key, mk in snp_table.items() if np.random.rand() <= keep_cds_probability and key in pb_cds_ids]
alpha_list = [dfe_alpha(cds_list, int(sample_size)) for sample_size in sample_size_list]
print(id_str)
print("Nbr of CDS={0}".format(len(cds_list)))
print("SFS sample size={0}".format(np.mean(sample_size_list)))
print("Proba to keep a CDS is {0}".format(keep_cds_probability))
x = [poly for alpha, omega, poly, div in alpha_list]
y = [div - omega for alpha, omega, poly, div in alpha_list]
plt.scatter(x, y)
for alpha, omega, poly, div in alpha_list:
    print("alpha={0}".format(alpha))
for alpha, omega, poly, div in alpha_list:
    print("omega={0}".format(omega))

plt.xlabel(r'$pn/ps')
plt.ylabel(r'$dn/ds-w_A^{NG}$')
plt.tight_layout()
plt.savefig("{0}/alpha_dfem_{1}_{2}_{3}.svg".format(data_path, version, GRCh, id_str), format="svg")
