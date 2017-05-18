from collections import defaultdict
import numpy as np
from math import floor
from collections import Counter
import os
import scipy.stats as st


def load_snp_table(data_path, version, GRCh, split=-1):
    snp_table = {}
    mk_data = open('{0}/{1}_{2}_estimates_snp{3}.out'.format(data_path, version, GRCh, "" if split == -1 else "_"+str(split)), 'r')
    mk_header = mk_data.readline().replace('\n', '').split('\t')
    for line in mk_data:
        line_split = line.replace('\n', '').split('\t')
        cds_id = line_split[0].split("_")[0]
        snp_table[cds_id] = dict(zip(mk_header[2:-2], [int(i) for i in line_split[2:-2]]))
        snp_table[cds_id][mk_header[1]] = line_split[1]
        assert len(line_split[1:-2]) == len(mk_header[1:-2])
        for index in [-1, -2]:
            if line_split[index] != "":
                snp_table[cds_id][mk_header[index]] = [int(i) for i in line_split[index].split(";")]
            else:
                snp_table[cds_id][mk_header[index]] = []
    mk_data.close()
    return snp_table

columns = sorted(["siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega", "pos"])


def str_to_table(table, label):
    return eval(label, {f: table[f] for f in columns})


def load_pb_table(data_path):
    dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
    return np.loadtxt("{0}/79_GRCh38_estimates_pb.out".format(data_path), dtype=dtype, skiprows=1)


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


def dfe_alpha(group, sample_size, dfe_path="/home/thibault/Tools/dfem", name="default"):
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
        range_sfs = range(1, int(floor(sample_size // 2))+1)
        assert len(range_sfs) * 2 == sample_size
        sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]

    sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
    sfs_filename = "sfs_{0}.txt".format(name)
    sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
    sfs_file.write("Homo\n")
    sfs_file.write("\t".join(sfs_list) + "\n")
    sfs_file.close()

    out_filename = "output_{0}.csv".format(name)
    model = "GammaExpo"
    os.system("{0}/dfem  -in {0}/{1} -out {0}/{2} -model {3}".format(dfe_path, sfs_filename, out_filename, model))
    out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
    headline = out_file.readline().split(',')
    index_model = headline.index("model")

    for line in out_file:
        line_split = line.split(',')
        if line_split[index_model] == model:
            omega_a = float(line_split[headline.index("omegaA")])
            out_file.close()
            return omega_a


def group_snp(group, cut_off=0):
    div_denom = sum([mk["Ds"] for mk in group]) * sum([mk["SitesN"] for mk in group])
    if div_denom == 0:
        return float("inf")
    div = sum([mk["Dn"] for mk in group]) * sum([mk["SitesS"] for mk in group]) / div_denom

    if cut_off == 0:
        poly_denom = sum([mk["Ps"] for mk in group]) * sum([mk["SitesN"] for mk in group])
        if poly_denom == 0:
            return float("inf")
        poly = sum([mk["Pn"] for mk in group]) * sum([mk["SitesS"] for mk in group]) / poly_denom
        return div-poly
    else:
        nc = int(cut_off*5008)
        poly_denom = sum([len([p for p in mk["SFSs"] if p > nc]) for mk in group]) * sum([mk["SitesN"] for mk in group])
        if poly_denom == 0:
            return float("inf")
        poly = sum([len([p for p in mk["SFSn"] if p > nc]) for mk in group]) * sum([mk["SitesS"] for mk in group]) / poly_denom
        return div - poly


def load_outliers(data_path, pb_table, snp_table):
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
    return outlier_set, non_outliers_list


def sfs_non_syn_and_syn(group, sample_size):
    sfs_n, sfs_s = [], []
    nbr_alleles = 5008
    sample_size = nbr_alleles if nbr_alleles < sample_size else sample_size
    if sample_size % 2 == 1:
        sample_size -= 1
    replacement = False
    for mk in group:
        sfs_s.extend(mk["SFSs"])
        sfs_n.extend(mk["SFSn"])

    sfs_s_sample = sample_sfs(sfs_s, sample_size, nbr_alleles, replacement)
    sfs_n_sample = sample_sfs(sfs_n, sample_size, nbr_alleles, replacement)

    return sfs_s_sample, sfs_n_sample