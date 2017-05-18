import os
from cds_libraries import build_dict_transcripts
import pickle as pickle
data_path = "/mnt/sda1/AdaptaPop/data"

gtf_folder = "archives_gtf"
gtf_path = "{0}/{1}".format(data_path, gtf_folder)
data_frame = []
error_file = open("{0}/gtf_errors.out".format(data_path), 'w')
error_full_file = open("{0}/gtf_errors_full.out".format(data_path), 'w')

for file in sorted(os.listdir(gtf_path), key=lambda x: int(x.split('.')[-2])):
    if file.endswith(".gtf"):
        error_file.write("\n\n" + file + "\n")
        error_full_file.write("\n\n" + file + "\n")
        cds_total = 0
        errors_cds_nf = []
        errors_cds_length = []

        dict_transcripts, not_confirmed_cds, full_transcripts = build_dict_transcripts(gtf_path, file)

        for tr_id, cds in dict_transcripts.items():
            cds_total += 1
            if not_confirmed_cds.get(tr_id):
                errors_cds_nf.append(tr_id)
                continue

            cds_length = cds.seq_length()
            if cds_length % 3 != 0:
                errors_cds_length.append(tr_id)

        cds_errors = len(errors_cds_nf)+len(errors_cds_length)
        if cds_errors > 0:
            error_file.write("{0} errors out of {1} coding sequences ({2:.2f}%)".format(cds_errors, cds_total, cds_errors*100./cds_total))
            error_full_file.write("{0} errors out of {1} coding sequences ({2:.2f}%)".format(cds_errors, cds_total, cds_errors*100./cds_total))
            if len(errors_cds_nf) > 0:
                error_file.write("\n\n{0} coding region start or end could not be confirmed ({1:.2f}%):\n".format(len(errors_cds_nf), len(errors_cds_nf)*100./cds_total))
                error_full_file.write("\n\n{0} coding region start or end could not be confirmed ({1:.2f}%):\n".format(len(errors_cds_nf), len(errors_cds_nf)*100./cds_total))
                # error_full_file.write("\n".join(["".join(full_transcripts[tr]) for tr in errors_cds_nf]))
            if len(errors_cds_length) > 0:
                error_file.write("\n\n{0} computed sequence size from exons starts and ends is not a multiple of 3 ({1:.2f}%):\n".format(len(errors_cds_length), len(errors_cds_length)*100./cds_total))
                error_full_file.write("\n\n{0} computed sequence size from exons starts and ends is not a multiple of 3 ({1:.2f}%):\n".format(len(errors_cds_length), len(errors_cds_length)*100./cds_total))
                error_full_file.write("\n".join(["".join(full_transcripts[tr]) for tr in errors_cds_length]))

        data_frame.append([file, len(errors_cds_nf), len(errors_cds_length), cds_total])

pickle.dump(data_frame, open("{0}/gtf_errors.p".format(data_path), "wb"))
error_file.close()
error_full_file.close()
print("Job completed")
