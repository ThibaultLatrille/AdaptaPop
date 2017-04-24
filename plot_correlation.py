import os
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from collections import Counter


def str_to_table(table, label, columns):
    return eval(label, {f: table[f] for f in columns})


def compute_mk_dict():
    mk_dict = {}
    mk_data = open('{0}/79_mk_test.out'.format(data_path), 'r')
    mk_header = mk_data.readline().replace('\n', '').split('\t')
    for line in mk_data:
        line_split = line.replace('\n', '').split('\t')
        mk_dict[line_split[0]] = dict(zip(mk_header[1:-2], [int(i) for i in line_split[1:-2]]))
        assert len(line_split[1:-2]) == len(mk_header[1:-2])
        for index in [-1, -2]:
            if line_split[index] != "":
                mk_dict[line_split[0]][mk_header[index]] = [int(i) for i in line_split[index].split(";")]
            else:
                mk_dict[line_split[0]][mk_header[index]] = []
    mk_data.close()
    return mk_dict


def alpha(group, version=2):
    if version == 0:
        denom = sum([mk["Dn"] for mk in group if mk]) * sum([mk["Ps"] for mk in group if mk])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] for mk in group if mk]) * sum([mk["Pn"] for mk in group if mk])
        return 1 - numer / denom
    elif version == 1:
        denom = sum([mk["Dn"] * mk["Ps"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk and mk["Ds"] + mk["Ps"] != 0)])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] * mk["Pn"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk and mk["Ds"] + mk["Ps"] != 0)])
        return 1 - numer / denom
    elif version == 2:
        dn = np.mean([mk["Dn"] for mk in group if mk])
        if dn == 0:
            return float("inf")
        ds = np.mean([mk["Ds"] for mk in group if mk])
        p = np.mean([mk["Ps"] / (mk["Ps"] + 1) for mk in group if mk])
        return 1 - ds * p / dn
    elif version == 3 or version == 4:
        return dfe_alpha(group, version)


def poly_div(group):
    mutation_dico = {}
    for snp in ["Pn", "Ps", "Dn", "Ds"]:
        mutation_dico[snp] = sum([mk[snp] for mk in group if mk])
    return mutation_dico


def omega_mean(group, mk_dict):
    seq_tot = 0
    omega_weighted = 0
    for cds_id, ome in group:
        omega_weighted += mk_dict[cds_id]["SeqLength"] * ome
        seq_tot += mk_dict[cds_id]["SeqLength"]
    return omega_weighted / seq_tot


def dfe_alpha(group, version):
    sites_n, sites_s, dn, ds, sfs_n_len, sfs_s_len, = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = Counter(), Counter()
    for mk in group:
        sites_n += mk["SitesN"]
        sites_s += mk["SitesS"]
        dn += mk["Dn"]
        ds += mk["Ds"]
        sfs_n_len += mk["SFSnLength"]
        sfs_s_len += mk["SFSsLength"]
        sfs_n.update(mk["SFSn"])
        sfs_s.update(mk["SFSs"])

    nbr_alleles = 2178
    if version == 3:
        dfe_path = "/home/thibault/Tools/dfe-alpha-2.15"
        sfs_filename = "sfs.txt"
        sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
        sfs_file.write("1\n{0}\n".format(nbr_alleles))
        for sfs, nbr_zero in [(sfs_s, sites_s - sfs_s_len), (sfs_n, sites_n - sfs_n_len)]:
            sfs_list = [str(nbr_zero)] + [str(sfs[i]) for i in range(1, nbr_alleles + 1)]
            for j in range(int(nbr_alleles / 2), nbr_alleles + 1):
                assert sfs_list[j] == "0"
            assert len(sfs_list) == nbr_alleles + 1
            sfs_file.write(" ".join(sfs_list) + "\n")
        sfs_file.close()

        base_config_str = "data_path_1\t{0}/data/\n" \
                          "data_path_2\tdata-three-epoch/\n" \
                          "sfs_input_file\t{0}/{1}\n" \
                          "est_dfe_results_dir\t{0}/{2}/\n"
        result_folder_str = "results_dir"
        config_filename = "config-est_dfe"
        opt_str = "site_class\t2\nfold\t1\nepochs\t1\nmean_s_variable\t1\nmean_s\t-0.1\nbeta_variable\t1\nbeta\t0.5"
        sel_config = open("{0}/{1}.txt".format(dfe_path, config_filename), 'w')
        sel_config.write(base_config_str.format(dfe_path, sfs_filename, result_folder_str) + opt_str)
        sel_config.close()

        os.system("{0}/est_dfe -c {0}/{1}.txt".format(dfe_path, config_filename))

        divergence_filename = "divergence.txt"
        divergence_file = open("{0}/{1}".format(dfe_path, divergence_filename), 'w')
        divergence_file.write("1 {0} {1}\n".format(sites_s - sfs_s_len, ds))
        divergence_file.write("0 {0} {1}\n".format(sites_n - sfs_n_len, dn))
        divergence_file.close()

        out_filename = "est_alpha_omega.out"
        est_alpha_config_filename = "config-est_alpha_omega"
        est_alpha_config_str = "data_path_1\t{0}/data/\n" \
                               "divergence_file\t{0}/{2}\n" \
                               "est_alpha_omega_results_file\t{0}/{3}\n" \
                               "est_dfe_results_file\t{0}/{1}/est_dfe.out\n" \
                               "neut_egf_file\t{0}/{1}/neut_egf.out\n" \
                               "sel_egf_file\t{0}/{1}/sel_egf.out\n" \
                               "do_jukes_cantor\t1\n" \
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
    elif version == 4:
        dfe_path = "/home/thibault/Tools/dfem"
        sfs_list = ["all_genes", str(nbr_alleles)]
        for sfs, nbr_site in [(sfs_s, sites_s), (sfs_n, sites_n)]:
            sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range(1, nbr_alleles + 1)]
        sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
        sfs_filename = "sfs.csv"
        sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
        sfs_file.write("#folded\n")
        sfs_file.write("\t".join(sfs_list) + "\n")
        sfs_file.close()

        out_filename = "sortie.csv"
        os.system(
            "{0}/dfem_omega  -in {0}/{1} -out {0}/{2} -model GammaExpo".format(dfe_path, sfs_filename, out_filename))

        out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
        index = out_file.readline().split(',').index("alpha")
        alpha_list = []
        for line in out_file:
            alpha_str = line.split(',')[index]
            if alpha_str != "NA":
                alpha_list.append(float(alpha_str))
        out_file.close()
        return sum(alpha_list) / len(alpha_list)


def bin_data(nbr_bins, n_plot, n_row, n_col, x_label, y_label, cds_and_omega_list, mk_dict, alpha_version):
    bin_size = int(len(cds_and_omega_list) / nbr_bins)

    grouped_cds_and_omega = [cds_and_omega_list[i:i + bin_size] for i in range(0, len(cds_and_omega_list), bin_size)]
    grouped_omega = [omega_mean(group, mk_dict) for group in grouped_cds_and_omega[2:-1]]
    grouped_mk_dict = [[mk_dict[cds_id] for cds_id, _ in group if mk_dict.get(cds_id)] for group in
                       grouped_cds_and_omega]
    grouped_alpha = [alpha(group, version=alpha_version) for group in grouped_mk_dict[2:-1]]

    grouped_alpha, grouped_omega = zip(*[(i, j) for i, j in zip(grouped_alpha, grouped_omega) if i <= 10.])

    plt.subplot(n_row, n_col, n_plot)
    plt.scatter(grouped_omega, grouped_alpha, linewidth=3)
    model = sm.OLS(grouped_alpha, grouped_omega)
    results = model.fit()
    a = results.params[0]
    idf = np.linspace(min(min(grouped_omega), 0), max(grouped_omega), 30)
    plt.plot(idf, idf, 'r--', label=r"$y=x$")
    plt.plot(idf, a * idf, '--', label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(min(min(grouped_omega), 0), max(grouped_omega))
    plt.legend()
    plt.title(r'${}$ bins'.format(len(grouped_omega)))


def plot_correlation(bins_list, params, alpha_version):
    assert alpha_version in [0, 1, 2, 3, 4]
    alpha_version_dict = {0: ("Fay", r'$\alpha_{Fay}$'),
                          1: ("TG", r'$\alpha_{TG}$'),
                          2: ("SEW", r'$\alpha_{SEW}$'),
                          3: ("EST", r'$\alpha_{EST}$'),
                          4: ("DFEM", r'$\alpha_{DFEM}$')}
    my_dpi = 96
    plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
    txt_file = open("{0}/79_correlation.out".format(data_path), 'w')
    n_plot = 1
    for col, x_label in params:

        omega_col = str_to_table(omega_table, col, columns)
        cds_id_table = [i[2:-1] for i in omega_table["CdsId"]]

        cds_and_omega_list = [(cds_id, omega) for cds_id, omega in zip(cds_id_table, omega_col) if cds_id in mk_dict]
        cds_and_omega_list.sort(key=lambda x: x[1])

        txt_file.write("CDS polymorphism data could not be found for these CDS:\n")
        txt_file.write(" ".join([cds_id for cds_id in cds_id_table if cds_id not in mk_dict]) + "\n")

        for nbr_bins in bins_list:
            bin_data(nbr_bins, n_plot, len(params), len(bins_list), x_label, alpha_version_dict[alpha_version][1],
                     cds_and_omega_list, mk_dict, alpha_version)
            n_plot += 1

    txt_file.close()
    plt.tight_layout()
    plt.savefig("{0}/correlation_alpha_{1}.svg".format(data_path, alpha_version_dict[alpha_version][0]), format="svg")


# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

columns = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega", "pos"])
dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)
mk_dict = compute_mk_dict()

params_alpha = [
    ("(siteomega-predmutsel)/siteomega", "$\\left< \\omega - \\omega_0 \\right> / \\left< \\omega \\right>$"),
    ("(mutselfreeomega-1)/mutselfreeomega", "$ (\\omega^* - 1) / \\omega^* $"),
    ("pos", "$\\left< P(\\omega > 1 ) \\right>$")]
params_omega = [("pos", "$\\left< P(\\omega > 1 ) \\right>$"),
                ("siteomega-predmutsel", "$\\left< \\omega - \\omega_0 \\right>$"),
                ("predmutselfreeomega*(mutselfreeomega-1)", "$\\left< \\omega_0^*  \\right>(\\omega^* - 1)$")]
bins_list = [50, 25, 10]

for alpha_version in [0, 1, 2, 3]:
    plot_correlation(bins_list, params_alpha, alpha_version)
    print('Plot completed')

print('Test completed')
