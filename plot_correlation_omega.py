import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

folders = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
data = "/pandata/tlatrill/AdaptaPop/data"
# data = "./data"

dtype = np.dtype([("tr_id", 'str')]+[(name, 'float64') for name in folders])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data), dtype=dtype, skiprows=1)

params = [("globalomega", "$\\omega$"),
          ("siteomega", "$\\left< \\omega_i \\right>$"),
          ("predmutsel", "$\\left< f_{i,j}/f_{i,k} \\right>$"),
          ("predmutselfreeomega", "$\\left< f_{i,j}^*/f_{i,k}^* \\right>$"),
          ("mutselfreeomega", "$\\omega^*$"),
          ("predmutselfreeomega*mutselfreeomega", "$\\omega^*\\left< f_{i,j}^*/f_{i,k}^* \\right>$")]
my_dpi = 96
plt.figure(figsize=(2 * 1920 / my_dpi, 2 * 1440 / my_dpi), dpi=my_dpi)


def prod(in_list):
    out = np.ones(len(in_list[0]))
    for i in in_list[0:]:
        out *= i
    return out


def access_column(table, col):
    list_params = col.split('*')
    if len(list_params) > 1:
        return prod([table[param] for param in list_params])
    else:
        list_params = col.split('+')
        if len(list_params) > 1:
            return sum([table[param] for param in list_params])
        else:
            return table[col]

n_plot = 1
for x, x_label in params:
    for y, y_label in params:
        plt.subplot(len(params), len(params), n_plot)
        n_plot += 1
        table_x = access_column(omega_table, x)
        table_y = access_column(omega_table, y)
        model = sm.OLS(table_y, table_x)
        results = model.fit()
        a = results.params[0]
        plt.scatter(table_x, table_y, linewidth=3)
        idf = np.linspace(0, max(table_x), 30)
        plt.plot(idf, idf)
        plt.plot(idf, a * idf, 'r-', label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared_adj))
        plt.legend()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim((0, max(table_x)))
        plt.ylim((0, max(table_y)))

plt.savefig("correlation_omega.svg", format="svg")
# plt.show()
print('Plot completed')
