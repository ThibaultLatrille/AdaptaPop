import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm


def str_to_table(table, label, columns):
    return eval(label, {f: table[f] for f in columns})


def alpha(group, version=2):
    if version == 0:
        denom = sum([mk["Dn"] for mk in group if mk]) * sum([mk["Ps"] for mk in group if mk])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] for mk in group if mk]) * sum([mk["Pn"] for mk in group if mk])
        return 1 - numer / denom
    elif version == 1:
        denom = sum([mk["Dn"] * mk["Ps"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk and mk["Ds"] + mk["Ps"] != 0)])
        if denom == 0:
            return float("inf")
        numer = sum([mk["Ds"] * mk["Pn"] / (mk["Ds"] + mk["Ps"]) for mk in group if (mk and mk["Ds"] + mk["Ps"] != 0)])
        return 1 - numer / denom
    elif version == 2:
        dn = np.mean([mk["Dn"] for mk in group if mk])
        if dn == 0:
            return float("inf")
        ds = np.mean([mk["Ds"] for mk in group if mk])
        p = np.mean([mk["Ps"] / (mk["Ps"] + 1) for mk in group if mk])
        return 1 - ds * p / dn


def poly_div(group):
    mutation_dico = {}
    for snp in ["Pn", "Ps", "Dn", "Ds"]:
        mutation_dico[snp] = sum([mk[snp] for mk in group if mk])
    return mutation_dico


def omega_mean(group, mk_dict):
    seq_tot = 0
    omega_weighted = 0
    for cds_id, ome in group:
        omega_weighted += mk_dict[cds_id]["SeqLength"] * ome
        seq_tot += mk_dict[cds_id]["SeqLength"]
    return omega_weighted / seq_tot


def bin_data(nbr_bins, n_plot, col, x_label, omega_list, mk_dict):
    bin_size = int(len(omega_list) / nbr_bins)
    txt_file.write("\n\nGrouping into a total of {0} bins, with {1} CDS per bin".format(nbr_bins, bin_size))
    grouped_omega = [omega_list[i:i + bin_size] for i in range(0, len(omega_list), bin_size)]

    grouped_alpha = [[mk_dict.get(cds_id) for cds_id, _ in group] for group in grouped_omega]
    for bin, (group_omega, group_alpha) in enumerate(zip(grouped_omega, grouped_alpha)):
        txt_file.write("\nBin {} containing {} CDS:\n".format(bin, len(group_omega)))
        txt_file.write("CdsId\t{0}\tPn\tPs\tDn\tDs\talpha\n".format(col))
        txt_file.write(
            "\n".join(["\t".join([cds_id, str(omega)] + [str(mk_dico[i]) for i in ["Pn", "Ps", "Dn", "Ds"]] + ["-"])
                       for (cds_id, omega), mk_dico in zip(group_omega, group_alpha)]) + "\n")
        dico_poly = poly_div(group_alpha)
        txt_file.write("\t".join(["Total", str(np.mean([omega for _, omega in group_omega]))]
                                 + [str(dico_poly[i]) for i in ["Pn", "Ps", "Dn", "Ds"]]
                                 + [str(alpha(group_alpha))]) + "\n")

    list_alpha_agglo = [alpha(group) for group in grouped_alpha[2:-1]]
    list_omega = [omega_mean(group, mk_dict) for group in grouped_omega[2:-1]]
    list_alpha_agglo, list_omega = zip(*[(i, j) for i, j in zip(list_alpha_agglo, list_omega) if i <= 1.])

    plt.subplot(2, 4, n_plot)
    plt.scatter(list_omega, list_alpha_agglo, linewidth=3)
    model = sm.OLS(list_alpha_agglo, list_omega)
    results = model.fit()
    a = results.params[0]
    idf = np.linspace(min(min(list_omega), 0), max(list_omega), 30)
    plt.plot(idf, idf, 'r--', label=r"$y=x$")
    plt.plot(idf, a * idf, '--', label=r"$y={0:.3g}x$ ($r^2={1:.3g})$".format(float(a), results.rsquared))
    plt.xlabel(x_label)
    plt.ylabel(r'$\alpha$')
    plt.xlim(min(min(list_omega), 0), max(list_omega))
    plt.legend()
    plt.title(r'${}$ bins'.format(len(list_omega)))


params = [("(siteomega-predmutsel)/siteomega", "$\\left< \\omega - \\omega_0 \\right> / \\left< \\omega \\right>$"),
          ("(mutselfreeomega-1)/mutselfreeomega", "$ (\\omega^* - 1) / \\omega^* $")]

# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

columns = sorted(["globalomega", "siteomega", "mutsel", "mutselfreeomega", "predmutsel", "predmutselfreeomega"])
dtype = np.dtype([("CdsId", 'str', 32)] + [(name, 'float64', 1) for name in columns])
omega_table = np.loadtxt("{0}/79_omega_estimated.out".format(data_path), dtype=dtype, skiprows=1)

mk_dict = {}
mk_data = open('{0}/79_mk_test.out'.format(data_path), 'r')
mk_header = mk_data.readline().replace('\n', '').split('\t')
for line in mk_data:
    line_split = line.replace('\n', '').split('\t')
    mk_dict[line_split[0]] = dict(zip(mk_header[1:-2], [int(i) for i in line_split[1:-2]]))
mk_data.close()

my_dpi = 96
plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)

txt_file = open("{0}/79_correlation.out".format(data_path), 'w')
n_plot = 1
for col, label in params:

    omega_col = str_to_table(omega_table, col, columns)
    cds_id_table = [i[2:-1] for i in omega_table["CdsId"]]

    omega_list = [(cds_id, omega) for cds_id, omega in zip(cds_id_table, omega_col) if cds_id in mk_dict]
    omega_list.sort(key=lambda x: x[1])

    txt_file.write("CDS polymorphism data could not be found for these CDS:\n")
    txt_file.write(" ".join([cds_id for cds_id in cds_id_table if cds_id not in mk_dict]) + "\n")

    for nbr_bins in [len(omega_list), 50, 25, 10]:
        bin_data(nbr_bins, n_plot, col, label, omega_list, mk_dict)
        n_plot += 1

txt_file.close()
plt.tight_layout()
plt.savefig("{0}/correlation_alpha.svg".format(data_path), format="svg")
print('Test completed')
