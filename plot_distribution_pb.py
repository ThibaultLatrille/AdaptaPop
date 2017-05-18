import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from cds_libraries import params_pb, load_pb_table, str_to_table

data_path = "./data"
pb_table = load_pb_table(data_path)

for negative in [True, False]:
    my_dpi = 96
    plt.figure(figsize=(2 * 1920 / my_dpi, 2 * 1440 / my_dpi), dpi=my_dpi)
    n_plot = 1
    rc_size = int(np.sqrt(len(params_pb))) + 1
    num_bins = 50
    for x, x_label in params_pb:
        plt.subplot(rc_size, rc_size, n_plot)
        n_plot += 1
        table_x = [val for val in str_to_table(pb_table, x) if val > 0 or negative]
        n, bins, patches = plt.hist(table_x, num_bins, normed=1, facecolor='green', alpha=0.5, label=r"$n={0}$".format(len(table_x)))
        y = mlab.normpdf(bins, np.mean(table_x), np.std(table_x))
        plt.plot(bins, y, 'r--')
        plt.xlabel(x_label)
        plt.ylabel('Probability')
        plt.xlim((min(table_x), max(table_x)))
    plt.legend()
    plt.tight_layout()
    plt.savefig("{0}/distribution_pb_{1}{2}.svg".format(data_path, num_bins, "" if negative else "_only_pos"), format="svg")
    plt.savefig("{0}/distribution_pb_{1}{2}.png".format(data_path, num_bins, "" if negative else "_only_pos"), format="png")

# plt.show()
print('Plot completed')
