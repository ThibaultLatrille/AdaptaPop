import os
from Bio import SeqIO
from Bio.pairwise2 import align, format_alignment
from pipeline.cds_libraries import build_dict_transcripts, build_dict_snps, codontable, complement
import numpy as np
from scipy.special import binom
from scipy import integrate

data_path = "/home/thibault/AdaptaPop/data"
homo_folder = "om_79_cds_homo"
filtered_path = "{0}/botero_fasta_filtered".format(data_path)
homo_path = "{0}/{1}".format(data_path, homo_folder)

phase3 = True
version = "88"
GRCh = "GRCh38"
sampling_effort = 40

if phase3:
    e_number_all, e_number_x, e_number_y = 5008, 3775, 1233
else:
    e_number_all, e_number_x, e_number_y = 2178, 1654, 524

amino_acids = sorted([aa for aa in set(codontable.values()) if aa != "stop"])
aa_index = {aa: k for k, aa in enumerate(amino_acids)}
assert (len(aa_index) == 20)
assert (max(aa_index.values()) == 19)
assert (min(aa_index.values()) == 0)


def sfs(i, n, f):
    precision = 10
    eps = 1e-5
    sample_sfs = np.zeros(n + 1)

    x_array = np.linspace(0 + eps, 1 - eps, 2 ** precision + 1)
    res_array = 2 * np.array([1 - np.exp(-f * (1 - x)) for x in x_array]) / (1 - np.exp(-f))
    for a in range(1, n + 1):
        y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, n - a - 1)
        if a == n:
            y_array[-1] = 2 * f / (1 - np.exp(-f))
        sample_sfs[a] = binom(n, a) * integrate.simps(y_array, x_array)

    sample_sfs /= np.nansum(sample_sfs)

    return np.nansum([x for a, x in enumerate(sample_sfs) if a >= i])


def nuc_seq_pos(seq, nbr_nucs):
    _pos = 0
    tmp_nbr_nucs = 0
    assert len(seq) >= nbr_nucs
    assert len(seq.replace("-", "")) >= nbr_nucs
    while True:
        if _pos >= len(seq):
            print("Last position for {0} nucleotides in sequence ".format(nbr_nucs) + seq)
            break
        if seq[_pos] != "-":
            tmp_nbr_nucs += 1
        if tmp_nbr_nucs == nbr_nucs:
            break
        _pos += 1
    return _pos


def build_site_profiles():
    profiles = {}
    folder_path = "{0}/pb_cleaned_mutsel".format(data_path)
    for file in os.listdir(folder_path):
        if file.endswith(".siteprofiles"):
            chain_name = file.replace('.siteprofiles', '').replace('_filtered_NT', '')
            cds_name = chain_name[:chain_name.rfind("_")]

            profil_dict = {}
            file_profiles = open("/".join([folder_path, file]), "r")
            file_profiles.readline()
            for line in file_profiles:
                l_split = line.split("\t")
                profil = list(map(float, l_split[1:]))
                profil_dict[int(l_split[0])] = profil
                assert len(profil) == 20
            file_profiles.close()
            if np.nansum([np.nansum(c) for c in profil_dict.values()]) > 0:
                profiles[cds_name] = profil_dict
    return profiles


print("Extracting .sitesprofiles files")
site_profiles = build_site_profiles()
print("Extracted {0} CDS with a profil".format(len(site_profiles)))

print("Extracting .gtf file")
dict_transcripts, not_confirmed_cds, full_transcripts = build_dict_transcripts(data_path,
                                                                               'Homo_sapiens_{0}_{1}.gtf'.format(
                                                                                   version, GRCh))
print("Extracted {0} CDS from .gtf".format(len(dict_transcripts)))

print("Extracting .vcf file")
dict_snps = build_dict_snps(data_path, 'Homo_sapiens_{0}_{1}_polymorphism_in_cds.vcf'.format(version, GRCh))
print("Extracted {0} SNPs from .vcf file".format(sum([len(b) for b in dict_snps])))

report_file = open('{0}/{1}_{2}_recent_selection.out'.format(data_path, version, GRCh), 'w')
tsv_file = open('{0}/{1}_{2}_recent_selection.tsv'.format(data_path, version, GRCh), 'w')
tsv_file.write(
    "CdsId\tChromosome\tSeqLength\tNbrExons\tRefCodon\tAltCodon\tRefAa\tAltAa\tRefNuc\tAltNuc\tAncNuc\tMinorNuc\tMAC\t"
    "AlleleCount\tSamplingEffort\tSelectionCoeff\tProba\n")

errors_cds = {"no_gtf": "Coding region not in the .gtf file",
              "not_confirmed": "Coding region start or end could not be confirmed",
              "no_inference": "Coding regions inference of Mutation-Selection model not performed",
              "raw_file_not_found": "Raw fasta file not found",
              "raw_homo_not_found": "Raw fasta homo sequence not found",
              "filtered_file_not_found": "Filtered fasta file not found",
              "filtered_homo_not_found": "Filtered fasta homo sequence not found",
              "stop_codon": "Coding region contains a stop codon in the sequence",
              "gtf_not_modulo_3": "The computed sequence size from exons starts and ends is not a multiple of 3",
              "fasta_not_modulo_3": "The sequence size is not a multiple of 3",
              "gtf_not_fasta": "The computed sequence size from exons starts and ends doesn't match the sequence size",
              "notaligned": "The raw and filtered sequences can't aligned",
              "misaligned": "The raw and filtered sequences are not properly aligned",
              "ali_too_large": "The alignment is larger than expected",
              "filtered_not_in_raw": "Nucleotide are in the filtered sequence but not in raw",
              "siteprofiles_not_filtered": "The filtered sequence size is different to siteprofiles"}

errors_cds_list = {k: [] for k in errors_cds.keys()}

errors_snp = {"fasta_not_vcf": "SNPs retrieved from the fasta are not equal to the reference",
              "ref_stop": "SNPs have stop codon as reference amino-acid",
              "alt_stop": "SNPs have stop codon as alternate amino-acid",
              "codon": "SNPs have non-identified reference or alternate amino-acid",
              "low_freq": "SNPS are not from 1000G project",
              "no_freq": "SNPS have no alleles count available",
              "no_ancestral": "SNPS don't have ancestral state",
              "ancestral_not_snp": "The ancestral is neither the reference nor the alternate",
              "not_in_filtered": "SNPS is not in the filtered sequence",
              "diff": "SNPS in the raw and filtered sequences are different"}
errors_snp_list = {k: [] for k in errors_snp.keys()}

cds_list = open('{0}/cds_homo_pan.out'.format(data_path), 'r')
cds_list_header = cds_list.readline().replace('\n', '').split('\t')
count_dict = {"cds": 0, "snp": 0, "snp_filtered": 0}


def analysis(line):
    line_split = line.replace('\n', '').split('\t')
    cds_dict = dict(zip(cds_list_header, line_split))
    tr_id = cds_dict["HomoTr"]
    cds_id = cds_dict["HomoCDS"]
    file_name = cds_dict["Filename"]
    count_dict["cds"] += 1
    if file_name != "ENSG00000110427_KIAA1549L":
        return

    if tr_id not in dict_transcripts:
        errors_cds_list["no_gtf"].append(cds_id)
        return

    if tr_id in not_confirmed_cds:
        errors_cds_list["not_confirmed"].append(file_name)
        return

    if cds_id not in site_profiles:
        errors_cds_list["no_inference"].append(file_name)
        return

    cds = dict_transcripts[tr_id]
    cds_length = cds.seq_length()
    if cds_length % 3 != 0:
        errors_cds_list["gtf_not_modulo_3"].append(file_name)
        return

    if os.path.isfile("{0}/{1}_raw_NT.fasta".format(homo_path, file_name)):
        fasta_seqs = SeqIO.parse(open("{0}/{1}_raw_NT.fasta".format(homo_path, file_name), 'r'), 'fasta')
        seqs = [fasta for fasta in fasta_seqs if fasta.id == "Homo"]
        if len(seqs) == 0:
            errors_cds_list["raw_homo_not_found"].append(file_name)
            return
        homo_seq = seqs[0].seq
        homo_seq_ungap = homo_seq.ungap('-')
    else:
        errors_cds_list["raw_file_not_found"].append(file_name)
        return

    t_homo_seq_ungap = len(homo_seq_ungap.translate(to_stop=True))
    if ((len(homo_seq_ungap) // 3) - t_homo_seq_ungap) != 1:
        print("stop_codon")
        errors_cds_list["stop_codon"].append(file_name)
        return

    seq_length = len(homo_seq_ungap)

    if seq_length % 3 != 0:
        print("fasta_not_modulo_3")
        errors_cds_list["fasta_not_modulo_3"].append(file_name)
        return

    if seq_length - cds_length != 3:
        print("gtf_not_fasta")
        errors_cds_list["gtf_not_fasta"].append(file_name)
        return

    filtered_cds_path = "{0}/allSpp_{1}_macse_NT_codon_gapClean40_Hmm5.fasta".format(filtered_path, cds_id)
    if os.path.isfile(filtered_cds_path):
        fasta_seqs = SeqIO.parse(open(filtered_cds_path, 'r'), 'fasta')
        seqs = [fasta for fasta in fasta_seqs if fasta.id == "Homo"]
        if len(seqs) == 0:
            errors_cds_list["filtered_homo_not_found"].append(file_name)
            return
        filtered_seq = str(seqs[0].seq)
    else:
        errors_cds_list["filtered_file_not_found"].append(file_name)
        return

    print("\n" + file_name)
    dict_match = {}
    letters = "ATCG-"

    for a in homo_seq:
        assert a in letters
    for a in filtered_seq:
        assert a in letters

    for f1 in letters:
        for f2 in letters:
            if f1 == "-" and f2 != "-":
                dict_match[(f1, f2)] = -500
            elif f1 != "-" and f2 != "-" and f1 == f2:
                dict_match[(f1, f2)] = 2
            elif f1 != "-" and f2 != "-" and f1 != f2:
                dict_match[(f1, f2)] = -10
            else:
                dict_match[(f1, f2)] = 0

    alignments = align.globaldd(homo_seq, filtered_seq, dict_match, -1000, -1000, -3, -0.5)

    if len(alignments) == 0:
        print("Raw and filtered could not be aligned")
        errors_cds_list["notaligned"].append(file_name)
        return

    best_alignment = alignments[0]

    if len(best_alignment[0]) != len(homo_seq):
        print("Raw sequence has been enlarged")
        errors_cds_list["ali_too_large"].append(file_name)
        return

    for s, first in enumerate(best_alignment[0]):
        if first == "-" and best_alignment[1][s] != "-":
            errors_cds_list["filtered_not_in_raw"].append(file_name)
            print("Misalignment position {0}".format(s))
            report_file.write(cds_id + "\n")
            report_file.write("Misalignment position {0}\n".format(s))
            report_file.write(format_alignment(*best_alignment) + "\n")
            return

    if len(site_profiles[cds_id]) != len(filtered_seq) // 3:
        print("Siteprofiles has a different size of filtered sequence")
        errors_cds_list["siteprofiles_not_filtered"].append(file_name)
        return

    snp_pass_nbr = 0
    snp_diff_nbr = 0
    if tr_id in dict_snps:
        for snp_id, chromosome, pos, ref_nt, alt_nt, meta_data in dict_snps[tr_id]:
            count_dict["snp"] += 1
            ref_aa, alt_aa, ref_codon, alt_codon = cds.amino_acid(homo_seq_ungap, pos, ref_nt, alt_nt)
            if ref_aa == '!':
                errors_snp_list["fasta_not_vcf"].append(
                    (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_nt, alt_nt, "0"))
            elif ref_aa == '-' or alt_aa == '-':
                errors_snp_list["codon"].append(
                    (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, "0"))
            elif ref_aa == 'stop':
                errors_snp_list["ref_stop"].append(
                    (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, "0"))
            elif alt_aa == 'stop':
                errors_snp_list["alt_stop"].append(
                    (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, "0"))
            else:
                if (";MAC=" in meta_data) and (";MAF=" in meta_data) and (";MA=" in meta_data):
                    meta_mac = meta_data[meta_data.index(";MAC=") + 5:].split(";")[0]
                    meta_maf = meta_data[meta_data.index(";MAF=") + 5:].split(";")[0]
                    e_number = round(int(meta_mac) / float(meta_maf))
                    if (e_number == e_number_all and (chromosome not in ["X", "Y", "MT"])) \
                            or (e_number == e_number_x and chromosome == "X") \
                            or (e_number == e_number_y and chromosome == "Y"):
                        lhs_homo_len = cds.nt_position(pos) + 1
                        homo_pos = nuc_seq_pos(best_alignment[0], lhs_homo_len)
                        if best_alignment[1][homo_pos] == "-":
                            errors_snp_list["diff"].append((file_name, snp_id, chromosome, str(pos), ref_codon,
                                                            alt_codon, ref_aa, alt_aa, meta_mac))
                            continue
                        lhs_filterd_len = len(best_alignment[1][:(homo_pos + 1)].replace("-", ""))
                        nt_filtered_pos = nuc_seq_pos(filtered_seq, lhs_filterd_len)
                        if nt_filtered_pos >= len(filtered_seq):
                            print("WATTTT the position is too high")
                            continue
                        elif nt_filtered_pos < 0:
                            print("WATTTT the position is too low")
                            continue
                        snp_pass_nbr += 1
                        if (cds.strand == "+" and filtered_seq[nt_filtered_pos] != ref_nt) or (
                                cds.strand == "-" and filtered_seq[nt_filtered_pos] != complement[ref_nt]):
                            errors_snp_list["diff"].append(
                                (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon,
                                 ref_aa, alt_aa, meta_mac))
                            snp_diff_nbr += 1
                        else:
                            count_dict["snp_filtered"] += 1
                            site_pos = int(nt_filtered_pos / 3) + 1
                            if site_pos > len(site_profiles[cds_id]):
                                print("Position of the site is {0} but sequence size is {1}".format(site_pos, len(
                                    site_profiles[cds_id])))
                            else:
                                if alt_aa != ref_aa:

                                    if ";AA=" not in meta_data:
                                        errors_snp_list["no_ancestral"].append(
                                            (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa,
                                             alt_aa, meta_mac))
                                        continue

                                    anc_nt = meta_data[meta_data.index(";AA=") + 4:].split(";")[0]
                                    if anc_nt not in [ref_nt, alt_nt]:
                                        errors_snp_list["ancestral_not_snp"].append(
                                            (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa,
                                             alt_aa, meta_mac))
                                        continue

                                    minor_nt = meta_data[meta_data.index(";MA=") + 4:].split(";")[0]
                                    assert minor_nt in [ref_nt, alt_nt]

                                    allele_count = round(float(meta_maf) * sampling_effort)
                                    delta_f = np.log(site_profiles[cds_id][site_pos][aa_index[alt_aa]])
                                    delta_f -= np.log(site_profiles[cds_id][site_pos][aa_index[ref_aa]])

                                    if anc_nt == alt_nt:
                                        delta_f = -delta_f

                                    if minor_nt == anc_nt:
                                        allele_count = sampling_effort - allele_count

                                    proba = sfs(allele_count, sampling_effort, delta_f)
                                    tsv_line = "\t".join(map(str, [file_name, snp_id, chromosome, str(pos),
                                                                   ref_codon, alt_codon, ref_aa, alt_aa,
                                                                   ref_nt, alt_nt, anc_nt, minor_nt, meta_mac,
                                                                   allele_count, sampling_effort, delta_f, proba]))
                                    if proba < 0.1:
                                        print(tsv_line)
                                        tsv_file.write(tsv_line)
                    else:
                        errors_snp_list["low_freq"].append(
                            (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon,
                             ref_aa, alt_aa, "0"))
                else:
                    errors_snp_list["no_freq"].append(
                        (file_name, snp_id, chromosome, str(pos), ref_codon, alt_codon, ref_aa, alt_aa, "0"))

    if snp_pass_nbr == 0 and len(dict_snps[tr_id]) != 0:
        print("Out of {0} total SNPs, none of them passed".format(len(dict_snps[tr_id])))
    elif snp_diff_nbr == 0:
        print("All of {0} SNPs which passed are right".format(snp_pass_nbr))
    else:
        print("{0} SNP out of {1} which passed are wrong".format(snp_diff_nbr, snp_pass_nbr))


for line in cds_list:
    if analysis(line) == -1:
        break

tsv_file.close()
cds_list.close()

error_file = open('{0}/{1}_{2}_recent_selection_errors.out'.format(data_path, version, GRCh), 'w')
error_light_file = open('{0}/{1}_{2}_recent_selection_errors_light.out'.format(data_path, version, GRCh), 'w')

errors_cds_total_count = sum([len(l) for l in errors_cds_list.values()])
if errors_cds_total_count > 0:
    tmp = "{0} errors out of {1} coding sequences ({2:.3f}%)".format(errors_cds_total_count, count_dict["cds"],
                                                                     errors_cds_total_count * 100. / count_dict["cds"])
    error_file.write(tmp)
    error_light_file.write(tmp)
    for key, msg in errors_cds.items():
        nbr_errors = len(errors_cds_list[key])
        if nbr_errors > 0:
            tmp = "\n\n{0} ({1} cds):\n".format(key, nbr_errors)
            error_file.write(tmp)
            error_light_file.write(tmp)
            error_file.write(" ".join(errors_cds_list[key]))

errors_snp_total_count = sum([len(l) for l in errors_snp_list.values()])
if errors_snp_total_count > 0:
    tmp = "\n\n{0} errors out of {1} SNPs ({2:.3f}%)".format(errors_snp_total_count, count_dict["snp"],
                                                             errors_snp_total_count * 100. / count_dict["snp"])
    error_file.write(tmp)
    error_light_file.write(tmp)
    for key, msg in errors_snp.items():
        nbr_errors = len(errors_snp_list[key])
        if nbr_errors > 0:
            tmp = "\n\n{0} {1} ({2:.3f}%):\n".format(nbr_errors, msg, 100 * nbr_errors / count_dict["snp"])
            error_file.write(tmp)
            error_light_file.write(tmp)
            error_file.write("CdsId\tSNPId\tChr\tPosition\tRefCodon\tAltCodon\tRefAA\tAltAA\tNbrAlleles\n")
            error_file.write("\n".join(["\t".join(er) for er in errors_snp_list[key]]))

error_file.close()
error_light_file.close()

print(count_dict["snp_filtered"])
print(count_dict["snp"])
if count_dict["snp"] != 0:
    print(count_dict["snp_filtered"] * 100 / count_dict["snp"])
else:
    print(0)
# assert snp_filtered + snp_errors == snp_total
print("Job completed")
