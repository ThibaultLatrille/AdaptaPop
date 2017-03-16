import os
from Bio import SeqIO
from lxml import etree
from cds_libraries import build_dict_transcripts, build_dict_snps, codontable

# data = "/pandata/tlatrill/AdaptaPop/data"
data = "./data"

tr_id = ""

dict_transcripts, not_confirmed_cds, full_transcripts = build_dict_transcripts(data, 'Homo_sapiens_79_GRCh37.gtf')
dict_snps = build_dict_snps(data, 'Homo_sapiens_79_polymorphism_in_cds.vcf')

path = "om_79_cds_homo"
txt_file = open(data + '/' + '79_mk_test.out', 'w')
txt_file.write("TrId\tSeqLength\tNbrExons\tPtot\tPn\tPs\tDtot\tDn\tDs\n")
cds_total = 0
errors_cds_nf = []
errors_cds_length = []
errors_cds_ungap_length = []
errors_cds_unequal_length = []
snp_total = 0
errors_snp_ref = []
errors_snp_stop = []
errors_snp_alt_stop = []
errors_snp_codon = []
for file in os.listdir(data + "/" + path):
    if file.endswith(".xml"):
        homo_seq_ungap, homo_seq, pan_seq = "", "", ""
        cds_total += 1
        file_name = file[:-4]

        root = etree.parse(data + "/" + path + "/" + file_name + ".xml").getroot()
        for specie in root.find('CDS').findall("speciesCDS"):
            if specie.attrib['species'] == 'Homo':
                tr_id = specie.find('infoCDS').find('ensidTr').text
                break

        if not_confirmed_cds.get(tr_id):
            errors_cds_nf.append(file_name)
            continue

        cds = dict_transcripts[tr_id]
        cds_length = cds.seq_length()
        if cds_length % 3 != 0:
            errors_cds_length.append(file_name)
            continue

        fasta_seqs = SeqIO.parse(open(data + "/" + path + "/" + file_name + "_raw_NT.fasta", 'r'), 'fasta')
        for fasta in fasta_seqs:
            if fasta.id == "Homo":
                homo_seq = fasta.seq
                homo_seq_ungap = homo_seq.ungap('-')
            elif fasta.id == "Pan":
                pan_seq = fasta.seq

        ungap_length = len(homo_seq_ungap)
        if ungap_length % 3 != 0:
            errors_cds_ungap_length.append(file_name)
            continue

        if cds_length != ungap_length - 3:
            errors_cds_unequal_length.append(file_name)
            continue

        if dict_snps.get(tr_id):
            list_aa_poly = []
            for snp_id, chromosome, pos, ref_nt, alt_nt in dict_snps[tr_id]:
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
                    list_aa_poly.append((ref_aa, alt_aa))

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
        txt_file.write(file_name + "\t" + str(len(homo_seq_ungap)) + "\t" + str(len(cds.exons)) + "\t"
                       + str(len(list_aa_poly)) + "\t" + pn + "\t" + ps + "\t"
                       + str(len(list_aa_div)) + "\t" + dn + "\t" + ds + "\n")
txt_file.close()

error_file = open(data + '/' + '79_mk_test_errors.out', 'w')
cds_errors = len(errors_cds_nf)+len(errors_cds_length)+len(errors_cds_ungap_length)+len(errors_cds_unequal_length)
if cds_errors > 0:
    error_file.write(str(cds_errors) + " errors out of " + str(cds_total) + " coding sequences ("
                     + '%.3f' % (cds_errors*100./cds_total) + "%)")
    if len(errors_cds_nf) > 0:
        error_file.write("\n\nCoding region start or end could not be confirmed ("
                         + str(len(errors_cds_nf)) + " cds):\n")
        error_file.write(" ".join(errors_cds_nf))
    if len(errors_cds_length) > 0:
        error_file.write("\n\nThe computed sequence size from exons starts and ends is not a multiple of 3 ("
                         + str(len(errors_cds_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_length))
    if len(errors_cds_ungap_length) > 0:
        error_file.write("\n\nThe sequence size is not a multiple of 3 ("
                         + str(len(errors_cds_ungap_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_ungap_length))
    if len(errors_cds_unequal_length) > 0:
        error_file.write("\n\nThe computed sequence size from exons starts and ends doesn't match the sequence size ("
                         + str(len(errors_cds_unequal_length)) + " cds):\n")
        error_file.write(" ".join(errors_cds_unequal_length))

snp_errors = len(errors_snp_ref)+len(errors_snp_codon)+len(errors_snp_stop)+len(errors_snp_alt_stop)
if snp_errors > 0:
    error_file.write("\n\n" + str(snp_errors) + " errors out of " + str(snp_total) + " SNPs ("
                     + '%.3f' % (snp_errors*100./snp_total) + "%)")
    if len(errors_snp_ref) > 0:
        error_file.write("\n\nSNPs retrieved from the fasta are not equal to the reference ("
                         + str(len(errors_snp_ref)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tSeqNt\tRefNt\tAltNt\tRefCodon\tAltCodon\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_ref]))
    if len(errors_snp_codon) > 0:
        error_file.write("\n\nThe reference or alternate amino-acid can't be identified ("
                         + str(len(errors_snp_codon)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_codon]))
    if len(errors_snp_stop) > 0:
        error_file.write("\n\nThe reference amino-acid is a stop codon ("
                         + str(len(errors_snp_stop)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_stop]))
    if len(errors_snp_alt_stop) > 0:
        error_file.write("\n\nThe alternate amino-acid is a stop codon ("
                         + str(len(errors_snp_alt_stop)) + " SNPs):\n")
        error_file.write("TrId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\n")
        error_file.write("\n".join(["\t".join(er) for er in errors_snp_alt_stop]))
error_file.close()

print("Job completed")
