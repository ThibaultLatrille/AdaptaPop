import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np


def str_to_table(table, label):
    return eval(label, {f: table[f] for f in columns})


cds_folder = "om_79_cds_mammals_no_pan_marsu"
columns = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"
cds_path = "{0}/{1}".format(data_path, cds_folder)
dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)

params = [("globalomega", "$\\omega$", False),
          ("siteomega", "$\\left< \\omega \\right>$", False),
          ("predmutsel", "$\\left< \\omega_0 \\right>$", False),
          ("predmutselfreeomega", "$\\left< \\omega_0^* \\right>$", False),
          ("mutselfreeomega", "$\\left< \\omega_* \\right>$", False),
          ("siteomega/predmutsel", "$\\left< \\omega / \\omega_0 \\right>$", False),
          ("siteomega-predmutsel", "$\\left< \\omega - \\omega_0 \\right>$", True),
          ("predmutselfreeomega*(mutselfreeomega-1)",
           "$ \\omega^*_0 \\left< \\omega^* - 1 \\right>$", True)]

my_dpi = 96
plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
n_plot = 1
rc_size = int(np.sqrt(len(params))) + 1
for num_bins in [25]:
    n_plot = 1
    for x, x_label, x_inter in params:
        plt.subplot(rc_size, rc_size, n_plot)
        n_plot += 1
        table_x = str_to_table(omega_table, x)
        n, bins, patches = plt.hist(table_x, num_bins, normed=1, facecolor='green', alpha=0.5)
        y = mlab.normpdf(bins, np.mean(table_x), np.std(table_x))
        plt.plot(bins, y, 'r--')
        plt.xlabel(x_label)
        plt.ylabel('Probability')
    plt.tight_layout()
    plt.savefig("{0}/distribution_omega_{1}.svg".format(data_path, num_bins), format="svg")

# plt.show()
print('Plot completed')
