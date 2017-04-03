import matplotlib.pyplot as plt
import numpy as np

folders = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"

dtype = np.dtype([("tr_id", 'str')]+[(name, 'float64') for name in folders])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data), dtype=dtype, skiprows=1)

params = ["globalomega", "siteomega", "predmutsel", "predmutselfreeomega", "mutselfreeomega"]
my_dpi = 96
plt.figure(figsize=(2 * 1920 / my_dpi, 2 * 1440 / my_dpi), dpi=my_dpi)

n_plot = 1
for x in params:
    for y in params:
        plt.subplot(len(params), len(params), n_plot)
        n_plot += 1
        plt.scatter(omega_table[x], omega_table[y], linewidth=3)
        idf = np.linspace(0, max(omega_table[x]), 30)
        plt.plot(idf, idf)
        plt.xlabel(x)
        plt.ylabel(y)
        plt.xlim((0, max(omega_table[x])))
        plt.ylim((0, max(omega_table[y])))

plt.savefig("correlation_omega.svg", format="svg")
plt.show()
print('Plot completed')
