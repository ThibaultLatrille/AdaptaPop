import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as st
from cds_libraries import load_pb_table, str_to_table, columns


def p_value(v, mean, std):
    if mean > v:
        return st.norm.cdf((v - mean) / std)
    else:
        return st.norm.cdf((mean - v) / std)

RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"

data_path = "./data"
pb_table = load_pb_table(data_path)

my_dpi = 128
fontsize = 16
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 2, 1)
table_x = str_to_table(pb_table, "predmutsel")
table_y = str_to_table(pb_table, "siteomega")

idf = np.linspace(min((min(table_x), 0)), max(table_x), 30)
plt.ylabel("$\\omega$", fontsize=fontsize)
plt.xlabel("$\\omega_0$", fontsize=fontsize)
plt.xlim((min((min(table_x), 0)), max(table_x)+0.01))
plt.ylim((min((min(table_y), 0)), max(table_y)+0.01))
model = sm.OLS(table_y, sm.add_constant(table_x))
results = model.fit()
b, a = results.params[0:2]
plt.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
             float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"), color=BLUE)
alpha = 0.005
prstd, iv_l, iv_u = wls_prediction_std(results, alpha=alpha)

txt_file = open('{0}/outliers.out'.format(data_path), 'w')
tr_ids = [i for i, y in enumerate(table_y) if (y < iv_l[i] or y > iv_u[i])]
txt_file.write("x={0}; y={1}".format("siteomega", "predmutsel"))
txt_file.write("\nCdsId\tx_{0}\ty_{1}\ty_predicted\ty_lower_bound\ty_upper_bound\tp_value\t".format("siteomega", "predmutsel") +
               "\t".join(columns) + "\n")
txt_file.write("\n".join(["\t".join(
    [str(j) for j in ([pb_table["CdsId"][i][2:-1],
                       table_x[i],
                       table_y[i],
                       results.fittedvalues[i],
                       iv_l[i],
                       iv_u[i],
                       p_value(table_y[i], results.fittedvalues[i], prstd[i])] +
                      [pb_table[col][i] for col in columns])
     ]) for i in tr_ids]) + "\n\n")
txt_file.close()

outlier_set = set()
cds_file = open('{0}/outliers.out'.format(data_path), 'r')
cds_file.readline()
cds_file.readline()
for line in cds_file:
    if line != '\n':
        outlier_set.add(line.split("\t")[0])
cds_file.close()

x_non_out, y_non_out = zip(*[(x, y) for i, (x, y) in enumerate(zip(table_x, table_y)) if pb_table["CdsId"][i][2:-1] not in outlier_set])
plt.scatter(x_non_out, y_non_out, linewidth=3, label=r"${0}$ non outliers CDS".format(len(x_non_out)), c=GREEN)

sorted_x = sorted(table_x)
plt.plot(sorted_x, sorted(iv_u), 'r--', color=BLUE)
plt.plot(sorted_x, sorted(iv_l), 'r--', label="{0}% confidence interval".format(100-alpha*100), color=BLUE)

x_out, y_out = zip(*[(x, y) for i, (x, y) in enumerate(zip(table_x, table_y)) if pb_table["CdsId"][i][2:-1] in outlier_set])
plt.scatter(x_out, y_out, linewidth=3, label=r"${0}$ outliers CDS".format(len(x_out)), c=RED)

plt.legend(fontsize=fontsize_legend)

plt.subplot(1, 2, 2)
table_omega_a = str_to_table(pb_table, "siteomega-predmutsel")
x_list, y_list = zip(*[(x, y) for i, (x, y) in enumerate(zip(table_x, table_omega_a)) if pb_table["CdsId"][i][2:-1] in outlier_set])
plt.scatter(x_list, y_list, linewidth=3, label=r"${0}$ outliers CDS".format(len(x_list)), c=RED)
idf = np.linspace(min((min(x_list), 0)), max(x_list), 30)
plt.xlabel("$\\omega_0$", fontsize=fontsize)
plt.ylabel("$\\omega_A =  \\omega - \\omega_0 $", fontsize=fontsize)
plt.xlim((min((min(x_list), 0)), max(x_list)+0.01))
plt.ylim((min(y_list)-0.01, max(y_list)+0.01))
model = sm.OLS(y_list, sm.add_constant(x_list))
results = model.fit()
b, a = results.params[0:2]
plt.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
             float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"), color=BLUE)
plt.legend(fontsize=fontsize_legend)
plt.tight_layout()
plt.savefig("{0}/figure_1.svg".format(data_path), format="svg")
plt.savefig("{0}/figure_1.png".format(data_path), format="png")

print('Plot completed')
