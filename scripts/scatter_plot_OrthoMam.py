import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import pandas as pd
import os

RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"
site = True
cds_list, table_omega_0, table_omega = [], [], []
for file in os.listdir("Experiments"):
    sitemutsel_path = "Experiments/" + file + "/sitemutsel_1.run.omegappgt1.000000"
    siteomega_path = "Experiments/" + file + "/siteomega_1.run.omegappgt1.000000"
    if not os.path.isfile(siteomega_path) or not os.path.isfile(sitemutsel_path):
        continue

    site_omega_0 = pd.read_csv(sitemutsel_path, sep="\t")["omega_0"].values
    site_omega = pd.read_csv(siteomega_path, sep="\t")["omega"].values
    if site:
        assert len(site_omega_0) == len(site_omega)
        table_omega_0.extend(site_omega_0)
        table_omega.extend(site_omega)
        cds_list.extend([file + str(i + 1) for i in range(len(site_omega))])
    else:
        table_omega_0.append(np.mean(site_omega_0))
        table_omega.append(np.mean(site_omega))
        cds_list.append(file)

table_omega_0 = np.array(table_omega_0)
table_omega = np.array(table_omega)
my_dpi = 128
fontsize = 16
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)
idf = np.linspace(min((min(table_omega_0), 0)), max(table_omega_0), 30)
plt.ylabel("$\\omega$", fontsize=fontsize)
plt.xlabel("$\\omega_0$", fontsize=fontsize)
xmin, xmax = min((min(table_omega_0), 0)), max(table_omega_0) + 0.01
ymin, ymax = min((min(table_omega), 0)), min(max(table_omega) + 0.01, 2.0)
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
pct = len([omega for omega_0, omega in zip(table_omega_0, table_omega) if omega > omega_0])
print("{0:3f}% have ω > ω0".format(100 * pct / len(table_omega)))

if site:
    # plt.hist2d(table_omega_0, table_omega, bins=100, range=[[0.0, 1.0], [0.0, 1.5]], norm=mpl.colors.LogNorm())
    plt.scatter(table_omega_0, table_omega)
    plt.plot(idf, idf)
else:
    model = sm.OLS(table_omega, sm.add_constant(table_omega_0))
    results = model.fit()
    b, a = results.params[0:2]
    plt.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
        float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"), color=BLUE)
    alpha = 0.005
    prstd, iv_l, iv_u = wls_prediction_std(results, alpha=alpha)

    cds_adaptive_list = [k for i, k in enumerate(cds_list) if table_omega[i] > iv_u[i]]
    cds_epistasis_list = [k for i, k in enumerate(cds_list) if table_omega[i] < iv_l[i]]
    cds_outliers_list = cds_adaptive_list + cds_epistasis_list
    x_non_out, y_non_out = zip(
        *[(table_omega_0[i], table_omega[i]) for i, k in enumerate(cds_list) if k not in cds_outliers_list])
    plt.scatter(x_non_out, y_non_out, linewidth=3, label=r"${0}$ non outliers CDS".format(len(x_non_out)), c=GREEN)

    sorted_x = sorted(table_omega_0)
    plt.plot(sorted_x, sorted(iv_u), 'r--', color=BLUE)
    plt.plot(sorted_x, sorted(iv_l), 'r--', label="{0}% confidence interval".format(100 - alpha * 100), color=BLUE)

    x_out, y_out = zip(*[(table_omega_0[i], table_omega[i]) for i, k in enumerate(cds_list) if k in cds_outliers_list])
    plt.scatter(x_out, y_out, linewidth=3, label=r"${0}$ outliers CDS".format(len(x_out)), c=RED)

    plt.legend(fontsize=fontsize_legend)

plt.tight_layout()
plt.savefig("scatterplot_OrthoMam{0}.pdf".format("_site" if site else ""), format="pdf")
plt.savefig("scatterplot_OrthoMam{0}.png".format("_site" if site else ""), format="png")

print('Plot completed')
