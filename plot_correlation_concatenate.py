import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

pbmpi_folder = "pb_output"
pbmpi_path = "{0}/{1}".format(data_path, pbmpi_folder)
concatenate_folder = "pb_concatenate"
concatenate_path = "{0}/{1}".format(data_path, concatenate_folder)

cds_omega = {}
bins = [10, 50]
cds_bins_omega = {i_bin: {} for i_bin in bins}
dtype = np.dtype([('#iter', 'int'), ('omega', 'float')])
cols = [0, 6]

for file in os.listdir(pbmpi_path):
    if file.endswith(".trace"):
        table = np.loadtxt("{0}/{1}".format(pbmpi_path, file), dtype=dtype, skiprows=1, usecols=cols)
        chain_name = file.replace('.trace', '').replace('_filtered_NT', '')
        cds_name = chain_name[:chain_name.rfind("_")]
        if cds_name not in cds_omega:
            cds_omega[cds_name] = []
        cds_omega[cds_name].append(np.mean(table['omega'][table['#iter'] > 100]))

pb_mean_omega = [(tr_id, np.mean(trace_ar)) for tr_id, trace_ar in cds_omega.items()]
pb_mean_omega.sort(key=lambda x: x[1])


for nbr_bins in bins:
    bin_size = int(len(pb_mean_omega) / nbr_bins)
    grouped_omega = [pb_mean_omega[i:i + bin_size] for i in range(0, len(pb_mean_omega), bin_size)]
    x_axis = []
    y_axis = []
    for i, group in enumerate(grouped_omega[:-1]):
        cds_bins_omega[nbr_bins][i] = []
        for chain in ["1", "2"]:
            file_name = "bins_{0}_group_{1}_{2}.trace".format(nbr_bins, i, chain)
            table = np.loadtxt("{0}/{1}".format(concatenate_path, file_name), dtype=dtype, skiprows=1, usecols=cols)
            cds_bins_omega[nbr_bins][i].append(np.mean(table['omega'][table['#iter'] > 50]))
        x_axis.extend([np.mean(cds_bins_omega[nbr_bins][i])]*len(group))
        y_axis.extend([omega for _, omega in group])
    idf = np.linspace(0, max(x_axis), 30)
    model = sm.OLS(y_axis, x_axis)
    results = model.fit()
    a = results.params[0]
    my_dpi = 96
    plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
    plt.scatter(x_axis, y_axis, linewidth=3)
    plt.plot(idf, idf)
    plt.plot(idf, a * idf, 'r-', label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared))
    plt.legend()
    plt.xlabel('$\\omega_{bin}$')
    plt.xlabel('$\\omega$')
    plt.xlim((0, max(x_axis)))
    plt.ylim((0, max(y_axis)))
    plt.tight_layout()
    plt.savefig("correlation_concatenate_{0}.svg".format(nbr_bins), format="svg")
    idf = np.linspace(0, max(y_axis), 30)
    model = sm.OLS(x_axis, y_axis)
    results = model.fit()
    a = results.params[0]
    my_dpi = 96
    plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
    plt.scatter(y_axis, x_axis, linewidth=3)
    plt.plot(idf, idf)
    plt.plot(idf, a * idf, 'r-', label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared))
    plt.legend()
    plt.xlabel('$\\omega_{bin}$')
    plt.xlabel('$\\omega$')
    plt.xlim((0, max(y_axis)))
    plt.ylim((0, max(x_axis)))
    plt.tight_layout()
    plt.savefig("correlation_concatenate_{0}_inv.svg".format(nbr_bins), format="svg")

print('Plot completed')
