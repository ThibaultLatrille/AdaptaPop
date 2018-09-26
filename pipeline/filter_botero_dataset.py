import os
from Bio import SeqIO
from Bio import Phylo
from collections import Counter


cleaned_folder = "alignmentsMammals_macseHmmclean"
ali_folder = "alignmentsMammals_macseHmmclean_ali"
tree = "mammals_no_homi_marsu.tree"
pb_file = "79_GRCh38_estimates_pb.out"

data_path = "/home/thibault/AdaptaPop/data"
cleaned_path = "{0}/{1}".format(data_path, cleaned_folder)
ali_path = "{0}/{1}".format(data_path, ali_folder)

species = [s.name for s in Phylo.read("{0}/{1}".format(data_path, tree), "newick").get_terminals()]
species_size = []

pb_path = open("{0}/{1}".format(data_path, pb_file), 'r')
pb_path.readline()
hit = 0
total = 0

for line in pb_path:
    cds_name = line.split("\t")[0]
    file_path = "{0}/allSpp_{1}_macse_NT_codon_gapClean40_Hmm5.fasta".format(cleaned_path, cds_name)
    total += 1
    if os.path.isfile(file_path):
        hit += 1
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

        ali_file = open("{0}/{1}.ali".format(ali_path, cds_name), 'w')
        ali_file.write("{0} {1}\n".format(len(filtered_seq), len(filtered_seq[0][1])))
        ali_file.write("\n".join([" ".join(id_seq) for id_seq in filtered_seq]))
        ali_file.close()

print("{0} hits".format(hit))
print("{0} total".format(total))
print('Job completed')
