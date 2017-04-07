import os
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as st

ops = {"+": (lambda x, y: x + y),
       "-": (lambda x, y: x - y),
       "*": (lambda x, y: x * y),
       "/": (lambda x, y: x / y)}


def str_to_table(table, label):
    i = max(label.find(op) for op in ops.keys())
    if i > -1:
        return ops[label[i]](str_to_table(table, label[:i]), str_to_table(table, label[i + 1:]))
    elif label.isdigit():
        return int(label)
    else:
        return table[label]


def p_value(v, mean, std):
    if mean > v:
        return st.norm.cdf((v - mean) / std)
    else:
        return st.norm.cdf((mean - v) / std)

cds_folder = "om_79_cds_mammals_no_pan_marsu"
outlier_folder = "outliers"
folders = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"
cds_path = "{0}/{1}".format(data_path, cds_folder)
file_name = "raw"
outlier_path = "{0}/{1}_{2}".format(data_path, outlier_folder, file_name)
os.system("rm -rf {0} && mkdir {0}".format(outlier_path))
dtype = np.dtype([("tr_id", 'str', 32)] + [(name, 'float64', 1) for name in folders])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)

if file_name == "adaptative":
    params = [("mutselfreeomega", "$\\left< \\omega^* \\right>$", False),
              ("siteomega/predmutsel", "$\\left< \\omega / \\omega_0 \\right>$", False),
              ("siteomega-predmutsel", "$\\left< \\omega - \\omega_0 \\right>$", True),
              ("predmutselfreeomega*mutselfreeomega-predmutselfreeomega",
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
txt_file = open('{0}/79_omega_outliers_{1}.out'.format(data_path, file_name), 'w')
for x, x_label, x_inter in params:
    for y, y_label, y_inter in params:
        plt.subplot(len(params), len(params), n_plot)
        n_plot += 1
        table_x = str_to_table(omega_table, x)
        table_y = str_to_table(omega_table, y)
        plt.scatter(table_x, table_y, linewidth=3)
        idf = np.linspace(min(min(table_x), 0), max(table_x), 30)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim((min(min(table_x), 0), max(table_x)))
        plt.ylim((min(min(table_y), 0), max(table_y)))
        if x_inter or y_inter:
            model = sm.OLS(table_y, sm.add_constant(table_x))
            results = model.fit()
            b, a = results.params[0:2]
            plt.plot(idf, a * idf + b, 'r-',
                     label=r"$y={0:.3g}x + {1:.3g}$ ($r^2={2:.3g})$".format(float(a), float(b), results.rsquared_adj))
        else:
            model = sm.OLS(table_y, table_x)
            results = model.fit()
            a = results.params[0]
            plt.plot(idf, idf)
            plt.plot(idf, a * idf, 'r-',
                     label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared_adj))
        if x != y:
            prstd, iv_l, iv_u = wls_prediction_std(results)
            sorted_x = sorted(table_x)
            plt.plot(sorted_x, sorted(iv_u), 'r--')
            plt.plot(sorted_x, sorted(iv_l), 'r--')
            tr_ids = [i for i, y in enumerate(table_y) if (y < iv_l[i] or y > iv_u[i])]
            txt_file.write("x={0}; y={1}".format(x, y))
            txt_file.write("\ntr_id\tx_{0}\ty_{1}\ty_predicted\ty_lower_bound\ty_upper_bound\tp_value\t".format(x, y) +
                           "\t".join(folders) + "\n")
            txt_file.write("\n".join(["\t".join(
                [str(j) for j in ([omega_table["tr_id"][i][2:-1],
                                   table_x[i],
                                   table_y[i],
                                   results.fittedvalues[i],
                                   iv_l[i],
                                   iv_u[i],
                                   p_value(table_y[i], results.fittedvalues[i], prstd[i])] +
                                  [omega_table[fold][i] for fold in folders])
                 ]) for i in tr_ids]) + "\n\n")
            for tr_id in [omega_table["tr_id"][i][2:-1] for i in tr_ids]:
                os.system("cp {0}/{1}_*.fasta {2}".format(cds_path, tr_id, outlier_path))
                os.system("cp {0}/{1}_*.rootree {2}".format(cds_path, tr_id, outlier_path))
                os.system("cp {0}/{1}.xml {2}".format(cds_path, tr_id, outlier_path))
        plt.legend()

txt_file.close()
plt.savefig("{0}/correlation_omega_{1}.svg".format(data_path, file_name), format="svg")
# plt.show()
print('Plot completed')
