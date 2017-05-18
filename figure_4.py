import numpy as np
import matplotlib as mpl
from math import floor
mpl.use('Agg')
import matplotlib.pyplot as plt
from cds_libraries import load_pb_table, load_snp_table, dfe_alpha, group_snp, load_outliers, sfs_non_syn_and_syn
import pickle as pickle
data_path = "/mnt/sda1/AdaptaPop/data"
version = "88"
GRCh = "GRCh38"
model = "GammaExpo"

snp_table = load_snp_table(data_path, version, GRCh, split=-1)
pb_table = load_pb_table(data_path)
outlier_set, non_outliers_list = load_outliers(data_path, pb_table, snp_table)

RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"
my_dpi = 128
fontsize = 16
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 1200 / my_dpi), dpi=my_dpi)

plt.subplot(2, 2, 1)
mk_distrib_size = 200000
mk_non_outliers = []
cut_off = 0.15
for _ in range(mk_distrib_size):
    group = [snp_table[cds_id] for cds_id in np.random.choice(non_outliers_list, len(outlier_set), replace=False)]
    mk_non_outliers.append(group_snp(group, cut_off=cut_off))

mk_outliers = group_snp([snp_table[cds_id] for cds_id in outlier_set], cut_off=cut_off)
above_outlier = len([w for w in mk_non_outliers if w > mk_outliers])
p_value = above_outlier / mk_distrib_size

num_bins = 150
weights = np.ones_like(mk_non_outliers)/float(len(mk_non_outliers))
a, _, _ = plt.hist(mk_non_outliers, num_bins, normed=False, weights=weights, facecolor=GREEN, histtype='stepfilled', alpha=0.75, label=r"Empirical distribution")

plt.xlabel("$\\omega_A^{MK}$", fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
y_max = 1.1*max(a)
plt.xlim((max(min(mk_non_outliers), -0.5), max(mk_non_outliers)))
plt.ylim((0, y_max))
plt.plot((mk_outliers, mk_outliers), (0, y_max), label=r"Outliers ($p_{\mathrm{value}}="+"{0:.3g}$)".format(p_value), color=RED)
plt.legend(fontsize=fontsize_legend)
print('MK')
print('Outlier w={0}'.format(mk_outliers))
print('Mean w={0}'.format(np.mean(mk_non_outliers)))
print('Std w={0}'.format(np.std(mk_non_outliers)))
print('n={0}'.format(mk_distrib_size))
print('nP(w>W)={0}'.format(above_outlier))
print('P(w>W)={0}'.format(p_value))
print('Subplot 1 done !')


plt.subplot(2, 2, 2)
f = open("sample_size_dico.p", 'rb')
sample_size_dico = pickle.load(f)
f.close()

sample_size_list = sorted(sample_size_dico.keys())

mean_dfem = np.array([np.mean(sample_size_dico[k]) for k in sample_size_list])
std_dfem = np.array([1.96 * np.std(sample_size_dico[k]) / np.sqrt(len(sample_size_dico[k])) for k in sample_size_list])
plt.plot(sample_size_list, mean_dfem, color=BLUE)
plt.plot(sample_size_list, mean_dfem+std_dfem, color=BLUE, linestyle='--')
plt.plot(sample_size_list, mean_dfem-std_dfem, color=BLUE, linestyle='--')
plt.scatter(sample_size_list, mean_dfem, c=GREEN)

plt.xlabel("Sampling effort", fontsize=fontsize)
plt.ylabel("$\\omega_A^{DFEM}$", fontsize=fontsize)
plt.xscale('log')
plt.xlim((min(sample_size_list), max(sample_size_list)))
plt.ylim((0, 0.6))
print('Subplot 2 done !')

plt.subplot(2, 2, 3)
sample_size = 24
bar_width = 0.35
opacity = 1

range_sfs = np.arange(1, int(floor(sample_size // 2)) + 1)
assert len(range_sfs) * 2 == sample_size

sfs_s_sample, sfs_n_sample = sfs_non_syn_and_syn([snp_table[cds_id] for cds_id in non_outliers_list], sample_size)
sfs_list = np.array([sfs_s_sample[i] for i in range_sfs])
y_non_outlier_s = sfs_list * 100 / sum(sfs_list)
plt.bar(range_sfs, y_non_outlier_s, bar_width, alpha=opacity, color=BLUE)
sfs_list = np.array([sfs_n_sample[i] for i in range_sfs])
y_non_outlier_n = sfs_list * 100 / sum(sfs_list)
plt.bar(range_sfs + bar_width, y_non_outlier_n, bar_width, alpha=opacity, color=LIGHTGREEN)

sfs_s_sample, sfs_n_sample = sfs_non_syn_and_syn([snp_table[cds_id] for cds_id in outlier_set], sample_size)
sfs_list = np.array([sfs_s_sample[i] for i in range_sfs])
y_outlier_s = - sfs_list * 100 / sum(sfs_list)
plt.bar(range_sfs, y_outlier_s, bar_width, alpha=opacity, color=BLUE, label='Synonymous')
sfs_list = np.array([sfs_n_sample[i] for i in range_sfs])
y_outlier_n = - sfs_list * 100 / sum(sfs_list)
plt.bar(range_sfs + bar_width, y_outlier_n, bar_width, alpha=opacity, color=LIGHTGREEN, label='Non-synonymous')

plt.ylim((1.1*min(min(y_outlier_n), min(y_outlier_s)), 1.1*max(max(y_non_outlier_n), max(y_non_outlier_s))))
plt.plot((0, max(range_sfs)+1), (0, 0), color='black')
plt.xlim((0, max(range_sfs)+1))
plt.xlabel("Minor allele frequency", fontsize=fontsize)
plt.ylabel('Proportion of SNPs (%)', fontsize=fontsize)
plt.legend()
print('Subplot 3 done !')


plt.subplot(2, 2, 4)
dfem_outliers = dfe_alpha([snp_table[cds_id] for cds_id in outlier_set], 12)
f = open("dfem_non_outliers.p", 'rb')
dfem_non_outliers = pickle.load(f)
f.close()

above_outlier = len([w for w in dfem_non_outliers if w > dfem_outliers])
if len(dfem_non_outliers) == 0:
    p_value = 0
else:
    p_value = above_outlier / len(dfem_non_outliers)
    plt.xlim((min(dfem_non_outliers), max(dfem_non_outliers)))

num_bins = 100
weights = np.ones_like(dfem_non_outliers)/float(len(dfem_non_outliers))
a, bins, patches = plt.hist(dfem_non_outliers, num_bins, normed=False, weights=weights, facecolor=GREEN,
                            histtype='stepfilled', alpha=0.75, label=r"Empirical distribution")

plt.xlabel("$\\omega_A^{DFEM}$", fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
y_max = 1.25*max(a)
plt.ylim((0, y_max))
plt.plot((dfem_outliers, dfem_outliers), (0, y_max), label=r"Outliers ($p_{\mathrm{value}}="+"{0:.3g}$)".format(p_value), color=RED)
plt.legend(fontsize=fontsize_legend)

print('DFEM')
print('Outlier w={0}'.format(dfem_outliers))
print('Mean w={0}'.format(np.mean(dfem_non_outliers)))
print('Std w={0}'.format(np.std(dfem_non_outliers)))
print('n={0}'.format(len(dfem_non_outliers)))
print('nP(w>W)={0}'.format(above_outlier))
print('P(w>W)={0}'.format(p_value))
print('Subplot 4 done !')

plt.savefig("{0}/figure_5.svg".format(data_path, num_bins), format="svg")
plt.savefig("{0}/figure_5.png".format(data_path, num_bins), format="png")
