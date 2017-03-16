import os
from Bio import SeqIO
from Bio import Phylo

path = "om_79_cds_mammals_no_pan_marsu"
tree = "mammals_no_homi_marsu.tree"
data = "/pandata/tlatrill/AdaptaPop/data"
# data = "./data"

species = [s.name for s in Phylo.read(data + "/" + tree, "newick").get_terminals()]


for file in os.listdir(data + "/" + path):
    if file.endswith(".fasta"):
        file_name = file[:-6]

        fasta_seqs = SeqIO.parse(open(data + "/" + path + "/" + file_name + ".fasta", 'r'), 'fasta')
        ali_file = open(data + "/" + path + "/" + file_name + ".ali", 'w')

        ids_seqs = [(fasta.id, str(fasta.seq)) for fasta in fasta_seqs if (fasta.id in species)]

        ali_file.write(str(len(ids_seqs)) + " " + str(len(ids_seqs[0][1])) + "\n")
        ali_file.write("\n".join([" ".join(id_seq) for id_seq in ids_seqs]))
        ali_file.close()

print('Job completed')
