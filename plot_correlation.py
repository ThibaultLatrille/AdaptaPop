import os
import matplotlib.pyplot as plt
import numpy as np

mk_file = "79_mk_test.out"
pbmpi_path = "pb_output"
# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"

mk_dict = {}
pb_dict = {}

mk_data = open(data + "/" + mk_file, 'r')
mk_header = mk_data.readline().replace('\n', '').split('\t')
for line in mk_data:
    line_split = line.replace('\n', '').split('\t')
    mk_dict[line_split[0]] = dict(zip(mk_header[1:], [int(i) for i in line_split[1:]]))
mk_data.close()


def alpha(group):
    denom = sum([mk["Dn"] for mk in group if mk]) * sum([mk["Ps"] for mk in group if mk])
    if denom == 0:
        return 1
    numer = sum([mk["Ds"] for mk in group if mk]) * sum([mk["Pn"] for mk in group if mk])
    return 1 - numer / denom


def poly_div(group):
    mutation_dico = {}
    for snp in ["Pn", "Ps", "Dn", "Ds"]:
        mutation_dico[snp] = sum([mk[snp] for mk in group if mk])
    return mutation_dico

for file in os.listdir(data + "/" + pbmpi_path):
    if file.endswith(".trace"):
        pb_data = open(data + "/" + pbmpi_path + "/" + file, 'r')
        chain_name = file.replace('.trace', '').replace('_filtered_NT', '')
        cds_name = chain_name[:chain_name.rfind("_")]
        line = pb_data.readline()
        pb_header = line.replace('\n', '').split('\t')
        for line in pb_data:
            split_line = line.replace('\n', '').split("\t")
            if int(split_line[0]) > 100:
                if pb_dict.get(cds_name):
                    for key, value in zip(pb_header, split_line):
                        pb_dict[cds_name][key].append(float(value))
                else:
                    pb_dict[cds_name] = dict(zip(pb_header, [[float(i)] for i in split_line]))
        pb_data.close()

txt_file = open(data + '/' + '79_correlation.out', 'w')
txt_file.write("CDS polymorphism data could not be found for these CDS:\n")
txt_file.write(" ".join([tr_id for tr_id in pb_dict.keys() if not mk_dict.get(tr_id)])+"\n")
pb_mean_omega = [(tr_id, np.mean(trace_dict["omega"])) for tr_id, trace_dict in pb_dict.items() if mk_dict.get(tr_id)]
pb_mean_omega.sort(key=lambda x: x[1])


def bin_data(nbr_bins, n_plot):
    bin_size = int(len(pb_mean_omega) / nbr_bins)
    txt_file.write("\n\nGrouping into a total of {0} bins, with {1} CDS per bin".format(nbr_bins, bin_size))
    grouped_omega = [pb_mean_omega[i:i + bin_size] for i in range(0, len(pb_mean_omega), bin_size)]

    grouped_alpha = [[mk_dict.get(tr_id) for tr_id, _ in group] for group in grouped_omega]
    for bin, (group_omega, group_alpha) in enumerate(zip(grouped_omega, grouped_alpha)):
        txt_file.write("\nBin {} containing {} CDS:\n".format(bin, len(group_omega)))
        txt_file.write("TrId\tomega\tPn\tPs\tDn\tDs\talpha\n")
        txt_file.write("\n".join(["\t".join([tr_id, str(omega)]+[str(mk_dico[i]) for i in ["Pn", "Ps", "Dn", "Ds"]]+["-"])
                                  for (tr_id, omega), mk_dico in zip(group_omega, group_alpha)]) + "\n")
        dico_poly = poly_div(group_alpha)
        txt_file.write("\t".join(["Total", str(np.mean([omega for _, omega in group_omega]))]
                                 + [str(dico_poly[i]) for i in ["Pn", "Ps", "Dn", "Ds"]]
                                 + [str(alpha(group_alpha))]) + "\n")
    grouped_omega = [[omega for _, omega in group] for group in grouped_omega]

    list_alpha_agglo = [alpha(group) for group in grouped_alpha]
    list_omega = [(np.mean(group) - 1) / np.mean(group) for group in grouped_omega]

    plt.subplot(2, 2, n_plot)
    plt.scatter(list_alpha_agglo, list_omega, linewidth=3)
    plt.ylabel(r'$\frac{\omega_* -1}{\omega_*}$')
    plt.xlabel(r'$\alpha$')
    plt.title(r'$\alpha$ are computed for each bin ({} bins)'.format(nbr_bins))

my_dpi = 96
plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
bin_data(50, 1)
bin_data(25, 2)
bin_data(10, 3)
bin_data(5, 4)
txt_file.close()

plt.show()
print('Test completed')
