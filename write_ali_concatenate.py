import os
from Bio import SeqIO
from Bio import Phylo
import numpy as np

ali_folder = "om_79_cds_mammals_no_pan_marsu"
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"

pbmpi_folder = "pb_output"
concatenate_folder = "concatenate"
concatenate_path = "{0}/{1}".format(data_path, concatenate_folder)
qsub_folder = "concatenate_qsub"
qsub_path = "{0}/{1}".format(data_path, qsub_folder)

tree = "mammals_no_homi_marsu.tree"
species = [s.name for s in Phylo.read(data_path + "/" + tree, "newick").get_terminals()]

pb_dict = {}

for file in os.listdir(data_path + "/" + pbmpi_folder):
    if file.endswith(".trace"):
        pb_data = open(data_path + "/" + pbmpi_folder + "/" + file, 'r')
        chain_name = file.replace('.trace', '')
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

pb_mean_omega = [(tr_id, np.mean(trace_dict["omega"])) for tr_id, trace_dict in pb_dict.items()]
pb_mean_omega.sort(key=lambda x: x[1])

for nbr_bins in [10, 50]:
    bin_size = int(len(pb_mean_omega) / nbr_bins)
    grouped_omega = [pb_mean_omega[i:i + bin_size] for i in range(0, len(pb_mean_omega), bin_size)]
    for i, group in enumerate(grouped_omega):
        ali_file = open("{0}/bins_{1}_group_{2}.ali".format(concatenate_path, nbr_bins, i), 'w')
        fasta_dico = {}
        for file_name, omega in group:
            fasta_seqs = SeqIO.parse(open("{0}/{1}/{2}.fasta".format(data_path, ali_folder, file_name), 'r'), 'fasta')
            for fasta in fasta_seqs:
                if fasta.id in species:
                    seq = str(fasta.seq)
                    if fasta_dico.get(fasta.id):
                        fasta_dico[fasta.id] += seq
                    else:
                        fasta_dico[fasta.id] = seq
        len_seq_list = [len(seq) for seq in fasta_dico.values()]
        assert len(set(len_seq_list)) == 1, "fail"

        ali_file.write("{0} {1}\n".format(len(len_seq_list), len_seq_list[0]))
        ali_file.write("\n".join(["{0} {1}".format(name, seq) for name, seq in fasta_dico.items()]))
        ali_file.close()

nbr_cpu = 4
os.system("rm {0}/*".format(qsub_path))
for file in os.listdir(concatenate_path):
    if file.endswith(".ali"):
        file_name = file[:-4]
        for chain in ["1", "2"]:
            qsub = open("{0}/{1}_{2}.pbs".format(qsub_path, file_name, chain), 'w')
            qsub.writelines("#!/bin/bash\n")
            qsub.writelines("#\n")
            qsub.writelines("#PBS -q q1day\n")
            qsub.writelines("#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu))
            qsub.writelines("#PBS -o /pandata/tlatrill/out_err/out_{0}_{1}\n".format(file_name, chain))
            qsub.writelines("#PBS -e /pandata/tlatrill/out_err/err_{1}_{1}\n".format(file_name, chain))
            qsub.writelines("#PBS -j oe\n")
            qsub.writelines("#PBS -W umask=022\n")
            qsub.writelines("#PBS -r n\n")
            qsub.writelines("#PBS -r n\n")
            qsub.writelines("TMP=/tmp/tlatrill$RANDOM\n")
            qsub.writelines("export TMPDIR=$TMP\n")
            qsub.writelines("mkdir -p {0}/pb_concatenate\n".format(data_path))
            command = "mpirun -n {0} ~/pbmpi/data/pb_mpi -f -s -x 1 400 -mutsel -freeomega -dp".format(nbr_cpu)
            command += " -d {0}/{1}".format(concatenate_path, file)
            command += " -T {0}/{1}".format(data_path, tree)
            command += " {0}/pb_concatenate/{1}_{2}\n".format(data_path, file_name, chain)
            qsub.writelines(command)
            qsub.writelines("rm -rf $TMP\n")
            qsub.writelines("rm {0}/{1}_{2}.pbs\n".format(qsub_path, file_name, chain))
            qsub.close()

        
print('Write completed')
