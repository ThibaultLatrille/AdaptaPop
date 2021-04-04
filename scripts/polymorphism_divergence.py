from collections import defaultdict
import numpy as np
from math import floor
from collections import Counter
import scipy.stats as st
import os
from scipy.special import binom
from scipy import integrate
import argparse
from Bio.Phylo.PAML import yn00
from Bio import SeqIO


def dNdS(homo_seq_fasta, pan_seq_fasta, tmp_path="./tmp"):
    homo_pan_str = "{0}/Homo_Pan.fasta".format(tmp_path)
    SeqIO.write([homo_seq_fasta, pan_seq_fasta], homo_pan_str, "fasta")
    yn = yn00.Yn00()
    yn.out_file = "{0}/yn.out".format(tmp_path)
    yn.alignment = homo_pan_str
    res = yn.run()["Homo"]["Pan"]['YN00']
    dn = int(round(res['N'] * res['dN']))
    ds = int(round(res['S'] * res['dS']))
    return dn, ds


def build_site_profiles():
    profiles = {}
    folder_path = "{0}/pb_cleaned_mutsel".format(data_path)
    for file in os.listdir(folder_path):
        if file.endswith(".siteprofiles"):
            chain_name = file.replace('.siteprofiles', '').replace('_filtered_NT', '')
            cds_name = chain_name[:chain_name.rfind("_")]

            profil_dict = {}
            file_profiles = open("/".join([folder_path, file]), "r")
            file_profiles.readline()
            for line in file_profiles:
                l_split = line.split("\t")
                profil = list(map(float, l_split[1:]))
                profil_dict[int(l_split[0])] = profil
                assert len(profil) == 20
            file_profiles.close()
            if np.nansum([np.nansum(c) for c in profil_dict.values()]) > 0:
                profiles[cds_name] = profil_dict
    return profiles


def build_dict_snps(_data_path, file_name):
    vcf_file = open("{0}/{1}".format(_data_path, file_name), 'r')
    _dict_snps = {}
    for line in vcf_file:
        if line[0] != '#':
            split_line = line.split("\t")
            if len(split_line[3]) == 1 and len(split_line[4]) == 1 and split_line[4] != ".":
                tr_id = split_line[11]
                if tr_id not in _dict_snps:
                    _dict_snps[tr_id] = []
                _dict_snps[tr_id].append(
                    (split_line[2], split_line[0], int(split_line[1]), split_line[3], split_line[4], split_line[7]))
    vcf_file.close()
    return _dict_snps


def load_pb_table(data_path):
    columns = sorted(["siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega", "pos"])

    dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
    return np.loadtxt("{0}/79_GRCh38_estimates_pb.out".format(data_path), dtype=dtype, skiprows=1)


def load_snp_table(data_path, version, GRCh, split=-1):
    snp_table = {}
    mk_data = open(
        '{0}/{1}_{2}_estimates_snp{3}.out'.format(data_path, version, GRCh, "" if split == -1 else "_" + str(split)),
        'r')
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


def sample_sfs(sfs, n_sample, n, replacement, theoretical=True):
    assert n_sample <= n
    fold = int(floor(n_sample // 2))
    if theoretical:
        out = {i: 0. for i in range(n_sample + 1)}
        for snp in sfs:
            x = np.arange(0, min(snp, n_sample) + 1)
            pdf = st.hypergeom(n, snp, n_sample).pmf(x)
            for k in x:
                if k < fold + 1:
                    out[k] += pdf[k]
                else:
                    out[n_sample - k] += pdf[k]
        assert sum(list(out.values())[fold + 1:]) == 0
        return out
    else:
        if replacement:
            sampled_sfs = np.random.binomial(n_sample, [i / n for i in sfs])
        else:
            sampled_sfs = np.random.hypergeometric(sfs, [n - i for i in sfs], n_sample)
        return Counter([i if i < fold + 1 else n_sample - i for i in list(sampled_sfs) if i > 0])


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
        return div - poly
    else:
        nc = int(cut_off * 5008)
        poly_denom = sum([len([p for p in mk["SFSs"] if p > nc]) for mk in group]) * sum([mk["SitesN"] for mk in group])
        if poly_denom == 0:
            return float("inf")
        poly = sum([len([p for p in mk["SFSn"] if p > nc]) for mk in group]) * sum(
            [mk["SitesS"] for mk in group]) / poly_denom
        return div - poly


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


def sfs(i, n, f):
    precision = 10
    eps = 1e-5
    sample_sfs = np.zeros(n + 1)

    x_array = np.linspace(0 + eps, 1 - eps, 2 ** precision + 1)
    res_array = 2 * np.array([1 - np.exp(-f * (1 - x)) for x in x_array]) / (1 - np.exp(-f))
    for a in range(1, n + 1):
        y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, n - a - 1)
        if a == n:
            y_array[-1] = 2 * f / (1 - np.exp(-f))
        sample_sfs[a] = binom(n, a) * integrate.simps(y_array, x_array)

    sample_sfs /= np.nansum(sample_sfs)

    return np.nansum([x for a, x in enumerate(sample_sfs) if a >= i])


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
        range_sfs = range(1, int(floor(sample_size // 2)) + 1)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--name', required=True, type=int,
                        dest="n", metavar="<qsub id>",
                        help="Qsub id")
    args = parser.parse_args()
    qsub_name = "subsampling_{0}".format(args.n)
    qsub_id = args.n

    model = "GammaExpo"
    data_path = "/panhome/tlatrill/AdaptaPop/data"
    dfe_path = "/panhome/tlatrill/AdaptaPop/dfem"
    GRCh = "GRCh38"
    version = "88"

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

    nbr_qsubs = 200
    np.random.seed(123456789)
    groups = [[snp_table[cds_id] for cds_id in np.random.choice(non_outliers_list, len(outlier_set), replace=False)] for
              _ in range(nbr_qsubs)]
    group = groups[qsub_id]

    sample_size_dico = {}
    sample_size_list = [int(i) for i in np.logspace(np.log10(8), np.log10(512), 2)]
    for sample_size in sample_size_list:
        sfs_filename = "sfs_{0}_{1}.txt".format(qsub_name, sample_size)
        out_filename = "out_{0}_{1}.csv".format(qsub_name, sample_size)
        dfem(group, dfe_path, sample_size, sfs_filename)

        os.system("{0}/dfem  -in {0}/{1} -out {0}/{2} -model {3}\n".format(dfe_path, sfs_filename, out_filename, model))

        out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
        headline = out_file.readline().split(',')
        index_model = headline.index("model")

        for line in out_file:
            line_split = line.split(',')
            if line_split[index_model] == model:
                omega_a = float(line_split[headline.index("omegaA")])
                if sample_size not in sample_size_dico:
                    sample_size_dico[sample_size] = []
                sample_size_dico[sample_size].append(omega_a)
        out_file.close()
        os.remove("{0}/{1}".format(dfe_path, out_filename))
        os.remove("{0}/{1}".format(dfe_path, sfs_filename))
