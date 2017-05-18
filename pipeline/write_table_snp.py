import os
from Bio import SeqIO
from cds_libraries import build_dict_transcripts, build_dict_snps
from Bio.Phylo.PAML import yn00
import numpy as np

data_path = "/mnt/sda1/AdaptaPop/data"
tmp_path = "/home/thibault/tmp"
homo_folder = "om_79_cds_homo"
homo_path = "{0}/{1}".format(data_path, homo_folder)

phase3 = True
version = "88"
GRCh = "GRCh38"
split = 1
np.random.seed(seed=123456789)

dict_transcripts, not_confirmed_cds, full_transcripts = build_dict_transcripts(data_path, 'Homo_sapiens_{0}_{1}.gtf'.format(version, GRCh))
dict_snps = build_dict_snps(data_path, 'Homo_sapiens_{0}_{1}_polymorphism_in_cds.vcf'.format(version, GRCh))

if phase3:
    e_number_all, e_number_x, e_number_y = 5008, 3775, 1233
else:
    e_number_all, e_number_x, e_number_y = 2178, 1654, 524

if version != "79":
    file_cds_homo_88 = "Homo_sapiens_{0}_{1}_cds_all".format(version, GRCh)
    fasta_dict = {}
    for fasta in SeqIO.parse(open("{0}/{1}.fasta".format(data_path, file_cds_homo_88), 'r'), 'fasta'):
        fasta_dict[fasta.id.split(".")[0]] = fasta.seq[:-3]

dict_mac = {}

txt_file = open('{0}/{1}_{2}_estimates_snp_{3}.out'.format(data_path, version, GRCh, split), 'w')
txt_file.write("CdsId\tChromosome\tSeqLength\tNbrExons\tPtot\tPn\tPs\tDtot\tDn\tDs\tSitesN\tSitesS\tSFSnLength\tSFSsLength\tSFSn\tSFSs\n")
cds_total = 0
errors_cds_nf = []
errors_cds_not_in_gtf = []
errors_cds_stop = []
errors_cds_length = []
errors_cds_ungap_length = []
errors_cds_unequal_length = []
snp_total = 0
snp_filtered = 0
errors_snp_ref = []
errors_snp_stop = []
errors_snp_alt_stop = []
errors_snp_codon = []
errors_snp_low_freq = []
errors_snp_no_freq = []


cds_list = open('{0}/cds_homo_pan.out'.format(data_path), 'r')
cds_list_header = cds_list.readline().replace('\n', '').split('\t')
for line in cds_list:
    line_split = line.replace('\n', '').split('\t')
    cds_dict = dict(zip(cds_list_header, line_split))
    tr_id = cds_dict["HomoTr"]
    cds_id = cds_dict["HomoCDS"]
    file_name = cds_dict["Filename"]
    homo_seq_ungap = ""
    cds_total += 1

    if tr_id not in dict_transcripts:
        errors_cds_not_in_gtf.append(cds_id)
        continue

    if tr_id in not_confirmed_cds:
        errors_cds_nf.append(file_name)
        continue

    cds = dict_transcripts[tr_id]
    cds_length = cds.seq_length()
    if cds_length % 3 != 0:
        errors_cds_length.append(file_name)
        continue

    if os.path.isfile("{0}/{1}_raw_NT.fasta".format(homo_path, file_name)):
        fasta_seqs = SeqIO.parse(open("{0}/{1}_raw_NT.fasta".format(homo_path, file_name), 'r'), 'fasta')
        for fasta in fasta_seqs:
            if fasta.id == "Homo" or fasta.id == "Pan":
                stop_index = max([fasta.seq.rfind(codon) for codon in ["TAA", "TAG", "TGA"]])
                fasta.seq = fasta.seq[:stop_index] + fasta.seq[stop_index+3:]
                if fasta.id == "Homo":
                    homo_seq_fasta = fasta
                    homo_seq_ungap = fasta.seq.ungap('-')
                if fasta.id == "Pan":
                    pan_seq_fasta = fasta
                    pan_seq_ungap = fasta.seq.ungap('-')
    else:
        print("{0}/{1}_raw_NT.fasta was not found !".format(homo_path, file_name))

    if len(homo_seq_ungap.translate(to_stop=True)) != len(homo_seq_ungap) // 3:
        errors_cds_stop.append(file_name)
        continue

    if len(pan_seq_ungap.translate(to_stop=True)) != len(pan_seq_ungap) // 3:
        errors_cds_stop.append(file_name)
        continue

    if version == "79":
        homo_seq = homo_seq_ungap
    else:
        homo_seq = fasta_dict[tr_id]

    seq_length = len(homo_seq)
    if seq_length % 3 != 0:
        errors_cds_ungap_length.append(file_name)
        continue

    if cds_length != seq_length:
        errors_cds_unequal_length.append(file_name)
        continue

    dict_mac[tr_id] = {"NS": [], "S": []}
    list_aa_poly = []
    if tr_id in dict_snps:
        for snp_id, chromosome, pos, ref_nt, alt_nt, meta_data in dict_snps[tr_id]:
            if split == int(round(np.random.rand())):
                snp_total += 1
                ref_aa, alt_aa, ref_codon, alt_codon = cds.amino_acid(homo_seq, pos, ref_nt, alt_nt)
                if ref_aa == '!':
                    errors_snp_ref.append((file_name, snp_id, chromosome, str(pos), alt_aa, ref_nt, alt_nt, ref_codon, alt_codon))
                elif ref_aa == '-' or alt_aa == '-':
                    errors_snp_codon.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif ref_aa == 'stop':
                    errors_snp_stop.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif alt_aa == 'stop':
                    errors_snp_alt_stop.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                else:
                    if ("MAC" in meta_data) and ("MAF" in meta_data):
                        meta_mac = meta_data[meta_data.index("MAC")+4:].split(";")[0]
                        meta_maf = meta_data[meta_data.index("MAF") + 4:].split(";")[0]
                        e_number = round(int(meta_mac) / float(meta_maf))
                        if (e_number == e_number_all and (chromosome not in ["X", "Y", "MT"])) \
                                or (e_number == e_number_x and chromosome == "X") \
                                or (e_number == e_number_y and chromosome == "Y"):
                            dict_mac[tr_id]["S" if ref_aa == alt_aa else "NS"].append(meta_mac)
                            snp_filtered += 1
                            list_aa_poly.append((ref_aa, alt_aa))
                        else:
                            errors_snp_low_freq.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, str(e_number)))
                    else:
                        errors_snp_no_freq.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))

    pn = str(len([0 for i, j in list_aa_poly if i != j]))
    ps = str(len([0 for i, j in list_aa_poly if i == j]))

    homo_seq_gap = homo_seq_fasta.seq.tomutable()
    pan_seq_gap = pan_seq_fasta.seq.tomutable()
    n = min([len(homo_seq_gap), len(pan_seq_gap)])
    for j in range(n//3):
        if split == int(round(np.random.rand())):
            homo_seq_gap[j * 3:j * 3 + 3] = "---"
            pan_seq_gap[j * 3:j * 3 + 3] = "---"

    homo_seq_fasta.seq = homo_seq_gap
    pan_seq_fasta.seq = pan_seq_gap

    homo_pan_str = "{0}/Homo_Pan.fasta".format(tmp_path)
    SeqIO.write([homo_seq_fasta, pan_seq_fasta], homo_pan_str, "fasta")
    yn = yn00.Yn00()
    yn.out_file = "{0}/yn.out".format(tmp_path)
    yn.alignment = homo_pan_str
    res = yn.run()["Homo"]["Pan"]['YN00']
    dn = int(round(res['N']*res['dN']))
    ds = int(round(res['S']*res['dS']))
    txt_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n".format(
        file_name, cds.chromosome, seq_length, len(cds.exons), len(list_aa_poly),
        pn, ps, dn + ds, dn, ds,
        int(round(res['N'])), int(round(res['S'])),
        len(dict_mac[tr_id]["NS"]), len(dict_mac[tr_id]["S"]),
        ";".join(dict_mac[tr_id]["NS"]), ";".join(dict_mac[tr_id]["S"])))
txt_file.close()
cds_list.close()

error_file = open('{0}/{1}_{2}_estimates_snp_{3}_errors.out'.format(data_path, version, GRCh, split), 'w')
error_light_file = open('{0}/{1}_{2}_estimates_snp_{3}_errors_light.out'.format(data_path, version, GRCh, split), 'w')

cds_errors = len(errors_cds_not_in_gtf)+len(errors_cds_nf)+len(errors_cds_stop)+len(errors_cds_length)+\
             len(errors_cds_ungap_length)+len(errors_cds_unequal_length)
if cds_errors > 0:
    tmp = "{0} errors out of {1} coding sequences ({2:.3f}%)".format(cds_errors, cds_total, cds_errors*100./cds_total)
    error_file.write(tmp)
    error_light_file.write(tmp)
    if len(errors_cds_not_in_gtf) > 0:
        tmp = "\n\nCoding region not in the .gtf file ({0} cds):\n".format(len(errors_cds_not_in_gtf))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_not_in_gtf))
    if len(errors_cds_nf) > 0:
        tmp = "\n\nCoding region start or end could not be confirmed ({0} cds):\n".format(len(errors_cds_nf))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_nf))
    if len(errors_cds_stop) > 0:
        tmp = "\n\nCoding region contains a stop codon in the sequence ({0} cds):\n".format(len(errors_cds_stop))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_stop))
    if len(errors_cds_length) > 0:
        tmp = "\n\nThe computed sequence size from exons starts and ends is not a multiple of 3 ({0} cds):\n".format(len(errors_cds_length))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_length))
    if len(errors_cds_ungap_length) > 0:
        tmp = "\n\nThe sequence size is not a multiple of 3 ({0} cds):\n".format(len(errors_cds_ungap_length))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_ungap_length))
    if len(errors_cds_unequal_length) > 0:
        tmp = "\n\nThe computed sequence size from exons starts and ends doesn't match the sequence size ({0} cds):\n".format(len(errors_cds_unequal_length))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_unequal_length))

snp_errors = len(errors_snp_ref)+len(errors_snp_codon)+len(errors_snp_stop)+len(errors_snp_alt_stop)+len(errors_snp_low_freq)+len(errors_snp_no_freq)
if snp_errors > 0:
    tmp = "\n\n{0} errors out of {1} SNPs ({2:.3f}%)".format(snp_errors, snp_total, snp_errors*100./snp_total)
    error_file.write(tmp)
    error_light_file.write(tmp)
    if len(errors_snp_ref) > 0:
        tmp = "\n\n{0} SNPs retrieved from the fasta are not equal to the reference ({1:.3f}%):\n".format(len(errors_snp_ref), 100*len(errors_snp_ref)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tSeqNt\tRefNt\tAltNt\tRefCodon\tAltCodon\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_ref]))
    if len(errors_snp_codon) > 0:
        tmp = "\n\n{0} SNPs have non-identified reference or alternate amino-acid ({1:.3f}%):\n".format(len(errors_snp_codon), 100*len(errors_snp_codon)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_codon]))
    if len(errors_snp_stop) > 0:
        tmp = "\n\n{0} SNPs have stop codon as reference amino-acid ({1:.3f}%):\n".format(len(errors_snp_stop), 100*len(errors_snp_stop)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_stop]))
    if len(errors_snp_alt_stop) > 0:
        tmp = "\n\n{0} SNPs have stop codon as alternate amino-acid ({1:.3f}%):\n".format(len(errors_snp_alt_stop), 100*len(errors_snp_alt_stop)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_alt_stop]))
    if len(errors_snp_no_freq) > 0:
        tmp = "\n\n{0} SNPS have no alleles count available ({1:.3f}%):\n".format(len(errors_snp_no_freq), 100*len(errors_snp_no_freq)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_no_freq]))
    if len(errors_snp_low_freq) > 0:
        tmp = "\n\n{0} SNPS are not from 1000G project ({1:.3f}%):\n".format(len(errors_snp_low_freq), 100*len(errors_snp_low_freq)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\tNbrAlleles\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_low_freq]))
error_file.close()
error_light_file.close()

print(snp_filtered)
print(snp_total)
print(snp_filtered * 100 / snp_total)
assert snp_filtered + snp_errors == snp_total
print("Job completed")
