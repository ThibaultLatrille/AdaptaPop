import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from math import floor
import numpy as np
import statsmodels.api as sm
from collections import Counter
from cds_libraries import load_snp_table, load_pb_table, str_to_table

phase3 = True
version = "88"
GRCh = "GRCh38"


def group_snp(group):
    div_denom = sum([mk["Ds"] for mk in group]) * sum([mk["SitesN"] for mk in group])
    if div_denom == 0:
        return float("inf")
    div = sum([mk["Dn"] for mk in group]) * sum([mk["SitesS"] for mk in group]) / div_denom
    poly_denom = sum([mk["Ps"] for mk in group]) * sum([mk["SitesN"] for mk in group])
    if poly_denom == 0:
        return float("inf")
    poly = sum([mk["Pn"] for mk in group]) * sum([mk["SitesS"] for mk in group]) / poly_denom
    return div-poly


def pb_mean(group_cds_and_pb, snp_table):
    seq_tot = 0
    pb_weighted = 0
    for cds_id, ome in group_cds_and_pb:
        pb_weighted += snp_table[cds_id]["SeqLength"] * ome
        seq_tot += snp_table[cds_id]["SeqLength"]
    return pb_weighted / seq_tot


def sample_sfs(sfs, n_sample, n, replacement):
    assert n_sample <= n
    if replacement:
        sampled_sfs = np.random.binomial(n_sample, [i/n for i in sfs])
    else:
        sampled_sfs = np.random.hypergeometric(sfs, [n - i for i in sfs], n_sample)
    return [i if i < n_sample // 2 + 1 else n_sample // 2 - i for i in list(sampled_sfs) if i > 0]


def dfem_snp(group, sample_size):
    sites_n, sites_s, dn, ds, sfs_n_len, sfs_s_len, = 0, 0, 0, 0, 0, 0
    sfs_n, sfs_s = Counter(), Counter()
    nbr_alleles = 5008 if phase3 else 2178
    sample_size = nbr_alleles if nbr_alleles < sample_size else sample_size
    replacement = False
    for mk in group:
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
    dfe_path = "/home/thibault/Tools/dfem"
    sfs_list = ["Summed", str(sample_size)]
    for sfs, nbr_site in [(sfs_n, sites_n), (sfs_s, sites_s)]:
        range_sfs = range(1, int(floor(sample_size // 2))+1)
        assert len(range_sfs) * 2 == sample_size
        sfs_list += [str(nbr_site)] + [str(sfs[i]) for i in range_sfs]
    sfs_list += [str(sites_n), str(dn), str(sites_s), str(ds)]
    sfs_filename = "sfs.txt"
    sfs_file = open("{0}/{1}".format(dfe_path, sfs_filename), 'w')
    sfs_file.write("Homo\n")
    sfs_file.write("\t".join(sfs_list) + "\n")
    sfs_file.close()

    model = "GammaExpo"  # GammaZero
    out_filename = "sortie.csv"
    os.system(
        "{0}/dfem -in {0}/{1} -out {0}/{2} -model {3}".format(dfe_path, sfs_filename, out_filename, model))

    out_file = open("{0}/{1}".format(dfe_path, out_filename), 'r')
    headline = out_file.readline().split(',')
    index_model = headline.index("model")
    for line in out_file:
        line_split = line.split(',')
        if line_split[index_model] == model:
            omega_a = float(line_split[headline.index("omegaA")])
            alpha = float(line_split[headline.index("alpha")])
            g_mean = float(line_split[headline.index(model + ":negGmean")])
            g_shape = float(line_split[headline.index(model + ":negGshape")])
            out_file.close()
            return omega_a, alpha, g_mean, g_shape


def plot_correlation(bins_list, pb_param, snp_table, pb_table, sample_size):
    snp_params = {("MK", r'$\omega_A^{MK}$'),
                  ("NG", r'$\omega_A^{NG}$')}
    for negative in [False, True]:
        txt_file = open('{0}/{1}_{2}_correlation_bins{3}.out'.format(data_path, version, GRCh, "" if negative else "_only_pos"), 'w')
        my_dpi = 96
        plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
        n_plot = 1
        for nbr_bins in bins_list:
            txt_file.write("Sample size = {}\n".format(sample_size))
            txt_file.write("{0} bins\n".format(nbr_bins))

            pb_col = str_to_table(pb_table, pb_param[0])
            cds_id_table = [i[2:-1] for i in pb_table["CdsId"]]

            cds_and_pb_list = [(cds_id, omega) for cds_id, omega in zip(cds_id_table, pb_col) if
                               (negative or omega > 0) and (cds_id in snp_table) and (snp_table[cds_id]["Chromosome"] not in ["X", "Y", "MT"])]
            cds_and_pb_list.sort(key=lambda x: x[1])

            bin_size = int(len(cds_and_pb_list) / nbr_bins)
            grouped_cds_and_pb = [cds_and_pb_list[i:i + bin_size] for i in range(0, len(cds_and_pb_list), bin_size)]
            grouped_pb = [pb_mean(group_cds_and_omega, snp_table) for group_cds_and_omega in grouped_cds_and_pb]
            grouped_snp_table = [[snp_table[cds_id] for cds_id, _ in group] for group in grouped_cds_and_pb]

            grouped_length = [sum([mk["SeqLength"] for mk in group]) for group in grouped_snp_table]

            txt_output_dico = {i: {"SeqLength": length, "PBomegaA": pb_val} for i, (pb_val, length) in enumerate(zip(grouped_pb, grouped_length))}
            for i, group in enumerate(grouped_snp_table):
                for stat in ["Dn", "Ds", "Pn", "Ps", "SitesS", "SitesN"]:
                    txt_output_dico[i][stat] = sum([mk[stat] for mk in group])

            for snp_version, y_label in snp_params:
                txt_file.write("{0} {1}\n".format(snp_version, y_label))

                if snp_version == "NG":
                    grouped_ng = [dfem_snp(group, sample_size) for group in grouped_snp_table]
                    grouped_snp = [omega_a for omega_a, _, _, _ in grouped_ng]
                    for i, (omega_a, alpha, g_mean, g_shape) in enumerate(grouped_ng):
                        txt_output_dico[i]["NGomegaA"] = omega_a
                        txt_output_dico[i]["alpha"] = alpha
                        txt_output_dico[i]["GammaMean"] = g_mean
                        txt_output_dico[i]["GammaShape"] = g_shape
                else:
                    grouped_snp = [group_snp(group) for group in grouped_snp_table]
                    for i, omega_a in enumerate(grouped_snp):
                        txt_output_dico[i]["MKomegaA"] = omega_a

                grouped_snp, grouped_pb = zip(*[(i, j) for i, j in zip(grouped_snp, grouped_pb) if i <= 1.])

                plt.subplot(len(bins_list), len(snp_params), n_plot)
                plt.scatter(grouped_pb, grouped_snp, linewidth=3)
                model = sm.OLS(grouped_snp, sm.add_constant(grouped_pb))
                results = model.fit()
                b, a = results.params[0:2]
                idf = np.linspace(min(min(grouped_pb), 0), max(grouped_pb), 30)
                plt.plot(idf, idf, 'r--', label=r"$y=x$")
                plt.plot(idf, a * idf + b, 'r-',
                         label=r"$y={0:.3g}x + {1:.3g}$ ($r^2={2:.3g})$".format(float(a), float(b), results.rsquared))
                plt.xlabel(pb_param[1])
                plt.ylabel(y_label)
                plt.xlim(min(min(grouped_pb), 0), max(grouped_pb))
                plt.legend()
                plt.title(r'${}$ bins'.format(len(grouped_pb)))

                n_plot += 1

            header = ["SeqLength", "PBomegaA", "MKomegaA", "NGomegaA", "alpha", "GammaMean", "GammaShape", "Dn", "Ds", "Pn", "Ps", "SitesS", "SitesN"]
            txt_file.write("Bin\t" + "\t".join(header) + "\n")
            for key, val in txt_output_dico.items():
                txt_file.write("\t".join([str(key)] + ["{0:.3g}".format(val[col]) for col in header])+"\n")
            txt_file.write("\n")

        plt.tight_layout()
        txt_file.close()
        plt.savefig('{0}/{1}_{2}_correlation_snp_pb{3}.svg'.format(data_path, version, GRCh, "" if negative else "_only_pos"), format="svg")
        plt.savefig('{0}/{1}_{2}_correlation_snp_pb{3}.png'.format(data_path, version, GRCh, "" if negative else "_only_pos"), format="png")


data_path = "/mnt/sda1/AdaptaPop/data"
pb_table = load_pb_table(data_path)
snp_table = load_snp_table(data_path, version, GRCh)

pb_param = ("siteomega-predmutsel", "$\\left< \\omega - \\omega_0 \\right>$")
bins_list = [12, 24]
sample_size = 48

plot_correlation(bins_list, pb_param, snp_table, pb_table, sample_size)

print('Test completed')
