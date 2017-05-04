import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from cds_libraries import load_snp_table, load_pb_table

data_path = "/mnt/sda1/AdaptaPop/data"
version = "88"
GRCh = "GRCh38"

pb_table = load_pb_table(data_path)
pb_cds_ids = set(cds_id[2:-1] for cds_id in pb_table["CdsId"])
snp_table = load_snp_table(data_path, version, GRCh)

header = "SeqLength\tSitesN\tSitesS\tPtot\tPn\tPs\tDtot\tDn\tDs".split("\t")

my_dpi = 96
plt.figure(figsize=(2*1920 / my_dpi, 2*1440 / my_dpi), dpi=my_dpi)
n_plot = 1
rc_size = int(np.sqrt(len(header)))
num_bins = 50
for x_label in header:
    plt.subplot(rc_size, rc_size, n_plot)
    n_plot += 1
    table_x = [mk[x_label] for key, mk in snp_table.items() if key in pb_cds_ids]
    n, bins, patches = plt.hist(table_x, num_bins, normed=1, facecolor='green', alpha=0.5, label=r"$n={0}$".format(len(table_x)))
    y = mlab.normpdf(bins, np.mean(table_x), np.std(table_x))
    plt.plot(bins, y, 'r--')
    plt.xlabel(x_label)
    plt.ylabel('Probability')
plt.legend()
plt.tight_layout()
plt.savefig("{0}/{1}_{2}_distribution_cds_{3}.svg".format(data_path, version, GRCh, num_bins), format="svg")
plt.savefig("{0}/{1}_{2}_distribution_cds_{3}.png".format(data_path, version, GRCh, num_bins), format="png")

# plt.show()
print('Plot completed')
