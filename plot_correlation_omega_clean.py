import os
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as st


def str_to_table(table, label):
    return eval(label, {f: table[f] for f in folders})


def p_value(v, mean, std):
    if mean > v:
        return st.norm.cdf((v - mean) / std)
    else:
        return st.norm.cdf((mean - v) / std)

cds_folder = "om_79_cds_mammals_no_pan_marsu"
file_name = "adaptative"
folders = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"
cds_path = "{0}/{1}".format(data_path, cds_folder)
dtype = np.dtype([("tr_id", 'str', 32)] + [(name, 'float64', 1) for name in folders])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)

if file_name == "adaptative":
    params = [("mutselfreeomega", "$\\left< \\omega^* \\right>$", False),
              ("siteomega/predmutsel", "$\\left< \\omega / \\omega_0 \\right>$", False),
              ("siteomega-predmutsel", "$\\left< \\omega - \\omega_0 \\right>$", True),
              ("predmutselfreeomega*(mutselfreeomega-1)",
               "$\\left< \\omega^*_0 (\\omega^* - 1) \\right>$", True)]
else:
    params = [("globalomega", "$\\omega$", False),
              ("siteomega", "$\\left< \\omega \\right>$", False),
              ("predmutsel", "$\\left< \\omega_0 \\right>$", False),
              ("predmutselfreeomega", "$\\left< \\omega_0^* \\right>$", False),
              ("mutselfreeomega*predmutsel", "$\\left<  \\omega^* \\omega_0 \\right>$", False),
              ("mutselfreeomega*predmutselfreeomega", "$\\left<  \\omega^* \\omega_0^* \\right>$", False)]

my_dpi = 96
plt.figure(figsize=(2 * 1920 / my_dpi, 2 * 1440 / my_dpi), dpi=my_dpi)
n_plot = 1
for x, x_label, x_inter in params:
    for y, y_label, y_inter in params:
        plt.subplot(len(params), len(params), n_plot)
        n_plot += 1
        table_x = str_to_table(omega_table, x)
        table_y = str_to_table(omega_table, y)
        if x_inter or y_inter:
            model = sm.OLS(table_y, sm.add_constant(table_x))
            results = model.fit()
        else:
            model = sm.OLS(table_y, table_x)
            results = model.fit()
        if x != y:
            prstd, iv_l, iv_u = wls_prediction_std(results)
            tr_ids = [i for i, y in enumerate(table_y) if iv_l[i] < y < iv_u[i]]
            x_filtered = [table_x[i] for i in tr_ids]
            y_filtered = [table_y[i] for i in tr_ids]
            plt.scatter(x_filtered, y_filtered, linewidth=3)
            idf = np.linspace(min(min(x_filtered), 0), max(x_filtered), 30)
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.xlim((min(min(x_filtered), 0), max(x_filtered)))
            plt.ylim((min(min(y_filtered), 0), max(y_filtered)))
            if x_inter or y_inter:
                model = sm.OLS(y_filtered, sm.add_constant(x_filtered))
                results = model.fit()
                b, a = results.params[0:2]
                plt.plot(idf, a * idf + b, 'r-',
                         label=r"$y={0:.3g}x + {1:.3g}$ ($r^2={2:.3g})$".format(float(a), float(b),
                                                                                results.rsquared_adj))
            else:
                model = sm.OLS(y_filtered, x_filtered)
                results = model.fit()
                a = results.params[0]
                plt.plot(idf, idf)
                plt.plot(idf, a * idf, 'r-',
                         label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared_adj))
            plt.legend()

plt.savefig("{0}/correlation_omega_cleaned_{1}.svg".format(data_path, file_name), format="svg")
# plt.show()
print('Plot completed')
