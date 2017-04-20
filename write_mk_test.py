import os
from Bio import SeqIO
from lxml import etree
from cds_libraries import build_dict_transcripts, build_dict_snps, codontable

# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"
homo_folder = "om_79_cds_homo"
homo_path = "{0}/{1}".format(data_path, homo_folder)
tr_id = ""

cut_off = 0.01

dict_transcripts, not_confirmed_cds, full_transcripts = build_dict_transcripts(data_path, 'Homo_sapiens_79_GRCh37.gtf')
dict_snps = build_dict_snps(data_path, 'Homo_sapiens_79_polymorphism_in_cds.vcf')
dict_freq = {}

txt_file = open('{0}/79_mk_test.out'.format(data_path), 'w')
txt_file.write("CdsId\tSeqLength\tNbrExons\tPtot\tPn\tPs\tDtot\tDn\tDs\tSFSnbr\tSFSn\tSFSs\n")
cds_total = 0
errors_cds_nf = []
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

for file in os.listdir(homo_path):
    if file.endswith(".xml"):
        homo_seq_ungap, homo_seq, pan_seq = "", "", ""
        cds_total += 1
        file_name = file[:-4]

        root = etree.parse("{0}/{1}.xml".format(homo_path, file_name)).getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_id = specie.find('infoCDS').find('ensidTr').text
                break

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
                if fasta.id == "Homo":
                    homo_seq = fasta.seq
                    homo_seq_ungap = homo_seq.ungap('-')
                elif fasta.id == "Pan":
                    pan_seq = fasta.seq
        else:
            print("{0}/{1}_raw_NT.fasta was not found !".format(homo_path, file_name))

        ungap_length = len(homo_seq_ungap)
        if ungap_length % 3 != 0:
            errors_cds_ungap_length.append(file_name)
            continue

        if cds_length != ungap_length - 3:
            errors_cds_unequal_length.append(file_name)
            continue

        dict_freq[tr_id] = {"NS": [], "S": []}
        list_aa_poly = []
        if tr_id in dict_snps:
            for snp_id, chromosome, pos, ref_nt, alt_nt, meta_data in dict_snps[tr_id]:
                snp_total += 1
                ref_aa, alt_aa, ref_codon, alt_codon = cds.amino_acid(homo_seq_ungap, pos, ref_nt, alt_nt)
                if ref_aa == '!':
                    errors_snp_ref.append((file_name, snp_id, chromosome, str(pos), alt_aa, ref_nt, alt_nt, ref_codon, alt_codon))
                elif ref_aa == '-' or alt_aa == '-':
                    errors_snp_codon.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif ref_aa == 'stop':
                    errors_snp_stop.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                elif alt_aa == 'stop':
                    errors_snp_alt_stop.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
                else:
                    if ("MA" in meta_data) and ("MAF" in meta_data) and ("AA" in meta_data):
                        meta_ma = meta_data[meta_data.index("MA")+3]
                        meta_aa = meta_data[meta_data.index("AA")+3]
                        meta_maf = meta_data[meta_data.index("MAF")+4:].split(";")[0]
                        freq = str(1-float(meta_maf)) if meta_ma == meta_aa else meta_maf
                        if float(freq) > cut_off:
                            dict_freq[tr_id]["S" if ref_aa == alt_aa else "NS"].append(freq)
                            snp_filtered += 1
                            list_aa_poly.append((ref_aa, alt_aa))
                        else:
                            list_aa_poly.append((ref_aa, alt_aa))
                            errors_snp_low_freq.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, freq))
                    else:
                        errors_snp_no_freq.append((file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa))
        list_aa_div = []
        for pos, (homo_n, pan_n) in enumerate(zip(homo_seq, pan_seq)):
            if homo_n != pan_n and homo_n != '-' and pan_n != '-':
                homo_nt_pos = len(homo_seq[:pos].ungap('-'))
                homo_frame = homo_nt_pos % 3
                homo_codon = str(homo_seq_ungap[homo_nt_pos - homo_frame:homo_nt_pos + 3 - homo_frame])
                pan_nt_pos = len(pan_seq[:pos].ungap('-'))
                pan_frame = pan_nt_pos % 3
                pan_codon = str(pan_seq.ungap('-')[pan_nt_pos - pan_frame:pan_nt_pos + 3 - pan_frame])
                list_aa_div.append((codontable[homo_codon], codontable[pan_codon]))

        pn = str(len([0 for i, j in list_aa_poly if i != j]))
        ps = str(len([0 for i, j in list_aa_poly if i == j]))
        dn = str(len([0 for i, j in list_aa_div if i != j and i != '-' and j != '-']))
        ds = str(len([0 for i, j in list_aa_div if i == j and i != '-' and j != '-']))
        txt_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(
            file_name, len(homo_seq_ungap), len(cds.exons), len(list_aa_poly),
            pn, ps, len(list_aa_div), dn, ds,
            len(dict_freq[tr_id]["NS"]) + len(dict_freq[tr_id]["S"]),
            ";".join(dict_freq[tr_id]["NS"]), ";".join(dict_freq[tr_id]["S"])))
txt_file.close()

error_file = open('{0}/79_mk_test_errors.out'.format(data_path), 'w')
error_light_file = open('{0}/79_mk_test_errors_light.out'.format(data_path), 'w')

cds_errors = len(errors_cds_nf)+len(errors_cds_length)+len(errors_cds_ungap_length)+len(errors_cds_unequal_length)
if cds_errors > 0:
    tmp = "{0} errors out of {1} coding sequences ({2:.3f}%)".format(cds_errors, cds_total, cds_errors*100./cds_total)
    error_file.write(tmp)
    error_light_file.write(tmp)
    if len(errors_cds_nf) > 0:
        tmp = "\n\nCoding region start or end could not be confirmed ({0} cds):\n".format(len(errors_cds_nf))
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write(" ".join(errors_cds_nf))
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
        tmp = "\n\n{0} SNPS have no frequency available ({1:.3f}%):\n".format(len(errors_snp_no_freq), 100*len(errors_snp_no_freq)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_no_freq]))
    if len(errors_snp_low_freq) > 0:
        tmp = "\n\n{0} SNPS have the mutated allele at low frequency ({1:.3f}%):\n".format(len(errors_snp_low_freq), 100*len(errors_snp_low_freq)/snp_total)
        error_file.write(tmp)
        error_light_file.write(tmp)
        error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\tFrequency\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_low_freq]))
error_file.close()
error_light_file.close()

print(snp_filtered)
print(snp_total)
print(snp_filtered * 100 / snp_total)
assert snp_filtered + snp_errors == snp_total
print("Job completed")
