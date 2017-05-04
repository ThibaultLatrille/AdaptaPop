import os
from math import floor
import numpy as np
from collections import Counter

phase3 = True


def sample_sfs(sfs, n_sample, n, replacement, cut_off):
    sfs_low = [i for i in list(sfs) if i > cut_off]
    assert n_sample <= n
    if replacement:
        sampled_sfs = np.random.binomial(n_sample, [i/n for i in sfs_low])
    else:
        sampled_sfs = np.random.hypergeometric(sfs_low, [n - i for i in sfs_low], n_sample)
    return [i if i < n_sample // 2 + 1 else n_sample // 2 - i for i in list(sampled_sfs) if i > 0]


def dfe_alpha(group, cut_off, sample_size, alpha_version, i, n):
    sites_n, sites_s, dn, ds, sfs_n_len, sfs_s_len, = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = Counter(), Counter()
    nbr_alleles = 5008 if phase3 else 2178
    sample_size = nbr_alleles if nbr_alleles < sample_size else sample_size
    replacement = False
    for mk in group:
        if mk["Chromosome"] not in ["X", "Y", "MT"]:
            sites_s += mk["SitesS"]
            sites_n += mk["SitesN"]
            ds += mk["Ds"]
            dn += mk["Dn"]
            sfs_s_sample = sample_sfs(mk["SFSs"], sample_size, nbr_alleles, replacement, cut_off)
            sfs_n_sample = sample_sfs(mk["SFSn"], sample_size, nbr_alleles, replacement, cut_off)
            if cut_off == 0 and sample_size == nbr_alleles and not replacement:
                assert len(sfs_s_sample) == len(mk["SFSs"])
                assert sum(sfs_s_sample) == sum(mk["SFSs"])
                for i, j in zip(sfs_s_sample, mk["SFSs"]):
                    assert i == j
            sfs_s.update(sfs_s_sample)
            sfs_n.update(sfs_n_sample)
            sfs_s_len += len(sfs_s_sample)
            sfs_n_len += len(sfs_n_sample)
    if alpha_version == 3:
        dfe_path = "/home/thibault/Tools/dfe-alpha-2.15"
        sfs_filename = "sfs.txt"
        sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
        sfs_file.write("1\n{0}\n".format(sample_size))
        for sfs, nbr_zero in [(sfs_n, sites_n - sfs_n_len), (sfs_s, sites_s - sfs_s_len)]:
            sfs_list = [str(nbr_zero)] + [str(sfs[i]) for i in range(1, sample_size + 1)]
            for j in range(int(nbr_alleles / 2) + 1, sample_size + 1):
                assert sfs_list[j] == "0"
            assert len(sfs_list) == sample_size + 1
            sfs_file.write(" ".join(sfs_list) + "\n")
        sfs_file.close()

        base_config_str = "data_path_1\t{0}/data/\n" \
                          "data_path_2\tdata-three-epoch/\n" \
                          "sfs_input_file\t{0}/{1}\n" \
                          "est_dfe_results_dir\t{0}/{2}/\n"
        result_folder_str = "results_dir"
        config_filename = "config-est_dfe"
        epochs = 1
        opt_str = "site_class\t2\nfold\t1\nepochs\t"
        if epochs == 1:
            opt_str += "1\n"
        else:
            opt_str += "2\nsearch_n2\t1\nt2_variable\t1\nt2\t50\n"
        opt_str += "mean_s_variable\t1\nmean_s\t-0.1\nbeta_variable\t1\nbeta\t0.5"
        sel_config = open("{0}/{1}.txt".format(dfe_path, config_filename), 'w')
        sel_config.write(base_config_str.format(dfe_path, sfs_filename, result_folder_str) + opt_str)
        sel_config.close()

        os.system("{0}/est_dfe -c {0}/{1}.txt".format(dfe_path, config_filename))

        divergence_filename = "divergence.txt"
        divergence_file = open("{0}/{1}".format(dfe_path, divergence_filename), 'w')
        divergence_file.write("1 {0} {1}\n".format(sites_n - sfs_n_len, dn))
        divergence_file.write("0 {0} {1}\n".format(sites_s - sfs_s_len, ds))
        divergence_file.close()

        out_filename = "est_alpha_omega.out"
        est_alpha_config_filename = "config-est_alpha_omega"
        est_alpha_config_str = "data_path_1\t{0}/data/\n" \
                               "divergence_file\t{0}/{2}\n" \
                               "est_alpha_omega_results_file\t{0}/{3}\n" \
                               "est_dfe_results_file\t{0}/{1}/est_dfe.out\n" \
                               "neut_egf_file\t{0}/{1}/neut_egf.out\n" \
                               "sel_egf_file\t{0}/{1}/sel_egf.out\n" \
                               "do_jukes_cantor\t0\n" \
                               "remove_poly\t0".format(dfe_path, result_folder_str,
                                                       divergence_filename, out_filename)
        est_alpha_config = open("{0}/{1}.txt".format(dfe_path, est_alpha_config_filename), 'w')
        est_alpha_config.write(est_alpha_config_str)
        est_alpha_config.close()

        os.system("{0}/est_alpha_omega -c {0}/{1}.txt".format(dfe_path, est_alpha_config_filename))

        out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
        out_line = out_file.readline().replace('\n', '')
        out_file.close()
        param_name = "alpha"  # "omega_A" | "alpha"
        return float(out_line[out_line.find(param_name) + len(param_name) + 1:].split(" ")[0])
    elif alpha_version == 5:
        sfs_list = ["{0}_{1}".format(n, i), str(sample_size)]
        for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
            range_sfs = range(1, int(floor(sample_size // 2))+1)
            assert len(range_sfs) * 2 == sample_size
            sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]
        sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
        return np.random.rand()


def group_snp(group, cut_off, sample_size, alpha_version, i, n):
    if alpha_version == 0:
        denom = sum([mk["Dn"] for mk in group if mk]) * sum([mk["Ps"] for mk in group])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] for mk in group if mk]) * sum([mk["Pn"] for mk in group])
        return 1 - numer / denom
    elif alpha_version == 1:
        denom = sum([mk["Dn"] * mk["Ps"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk["Ds"] + mk["Ps"] != 0)])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] * mk["Pn"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk["Ds"] + mk["Ps"] != 0)])
        return 1 - numer / denom
    elif alpha_version == 2:
        dn = np.mean([mk["Dn"] for mk in group])
        if dn == 0:
            return float("inf")
        ds = np.mean([mk["Ds"] for mk in group])
        p = np.mean([mk["Ps"] / (mk["Ps"] + 1) for mk in group])
        return 1 - ds * p / dn
    else:
        return dfe_alpha(group, cut_off, sample_size, alpha_version, i=i, n=n)