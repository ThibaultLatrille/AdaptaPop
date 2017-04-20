import os
from Bio import SeqIO
from Bio import Phylo
from collections import Counter

cleaned_folder = "botero_dataset"
homo_folder = "om_79_cds_homo"
ali_folder = "botero_ali_filtered"
tree = "mammals_no_homi_marsu.tree"
data_path = "/pandata/tlatrill/AdaptaPop/data"
cleaned_path = "{0}/{1}".format(data_path, cleaned_folder)
homo_path = "{0}/{1}".format(data_path, homo_folder)
ali_path = "{0}/{1}".format(data_path, ali_folder)

os.system("rm {0}/qsub/*".format(data_path))

species = [s.name for s in Phylo.read("{0}/{1}".format(data_path, tree), "newick").get_terminals()]
species_size = []
discarded = 0
count = 0
nbr_cpu = 4
minimum_nbr_species = 16
for file in os.listdir(homo_path):
    if file.endswith(".fasta") and ('raw_NT' in file):
        cds_name = file.split("_")[0]
        file_path = "{0}/allSpp_{1}_macse_NT_codon_gapClean40_Hmm5.fasta".format(cleaned_path, cds_name)
        if os.path.isfile(file_path):
            fasta_seqs = SeqIO.parse(open(file_path, 'r'), 'fasta')

            ids_seqs = [(fasta.id, str(fasta.seq)) for fasta in fasta_seqs]

            len_set = set([len(seq) for id_seq, seq in ids_seqs])
            assert len(len_set) == 1, "Sequences have different length"
            len_seq = len_set.pop()
            assert len_seq % 3 == 0, "Sequences are not multiple of 3"
            species_size.append(len(ids_seqs))

            for counter, len_seq in [(Counter(seq), len(seq)) for id_seq, seq in ids_seqs]:
                assert counter['-'] != len_seq, "Empty sequence"

            filtered_seq = [(id_seq, seq) for id_seq, seq in ids_seqs if (id_seq in species)]
            if len(filtered_seq) > minimum_nbr_species:
                ali_file = open("{0}/{1}.ali".format(ali_path, cds_name), 'w')
                ali_file.write("{0} {1}\n".format(len(filtered_seq), len(filtered_seq[0][1])))
                ali_file.write("\n".join([" ".join(id_seq) for id_seq in filtered_seq]))
                ali_file.close()
                for option, opt_name in [("-siteomega", "siteomega"),
                                         ("-mutsel -dp", "mutsel"),
                                         ("-mutsel -freeomega -dp", "mutselfreeomega")]:
                    for chain in ["1"]:
                        count += 1
                        qsub = open("{0}/qsub/batch_{1}_{2}_{3}_{4}.pbs".format(data_path, count//500, cds_name, opt_name, chain), 'w')
                        qsub.writelines("#!/bin/bash\n")
                        qsub.writelines("#\n")
                        qsub.writelines("#PBS -q q1day\n")
                        qsub.writelines("#PBS -l nodes=1:ppn={0},mem=12gb\n".format(nbr_cpu))
                        qsub.writelines("#PBS -o /pandata/tlatrill/out_err/out_{0}_{1}_{2}\n".format(cds_name, opt_name, chain))
                        qsub.writelines("#PBS -e /pandata/tlatrill/out_err/err_{0}_{1}_{2}\n".format(cds_name, opt_name, chain))
                        qsub.writelines("#PBS -j oe\n")
                        qsub.writelines("#PBS -W umask=022\n")
                        qsub.writelines("#PBS -r n\n")
                        qsub.writelines("#PBS -r n\n")
                        qsub.writelines("TMP=/tmp/tlatrill$RANDOM\n")
                        qsub.writelines("export TMPDIR=$TMP\n")
                        qsub.writelines("mkdir -p {0}/pb_cleaned_{1}\n".format(data_path, opt_name))
                        command = "mpirun -n {0} ~/pbmpi2/data/pb_mpi -f -s -x 1 600 {1}".format(nbr_cpu, option)
                        command += " -d {0}/{1}.ali".format(ali_path, cds_name)
                        command += " -T {0}/{1}".format(data_path, tree)
                        command += " {0}/pb_cleaned_{1}/{2}_{3}\n".format(data_path, opt_name, cds_name, chain)
                        qsub.writelines(command)
                        qsub.writelines("rm -rf $TMP\n")
                        qsub.writelines("rm {0}/qsub/batch_{1}_{2}_{3}_{4}.pbs\n".format(data_path, count//500, cds_name, opt_name, chain))
                        qsub.close()
            else:
                discarded += 1
        else:
            discarded += 1

print(len(species_size))
print(Counter(species_size))
print(discarded)
print('Job completed')
