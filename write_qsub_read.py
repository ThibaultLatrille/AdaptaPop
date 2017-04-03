import os

ali_folder = "om_79_cds_mammals_no_pan_marsu_1"
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"

ali_path = "{0}/{1}".format(data_path, ali_folder)
os.system("rm {0}/qsub/*".format(data_path))
for file in os.listdir(ali_path):
    if file.endswith("filtered_NT.ali"):
        file_name = file[:-4]
        for chain in ["1", "2"]:
            qsub = open("{0}/qsub/{1}_{2}.pbs".format(data_path, file_name, chain), 'w')
            qsub.writelines("#!/bin/bash\n")
            qsub.writelines("#\n")
            qsub.writelines("#PBS -q q1day\n")
            qsub.writelines("#PBS -l nodes=1:ppn=1,mem=1gb\n")
            qsub.writelines("#PBS -o /pandata/tlatrill/read/read_out_{0}_{1}\n".format(file_name, chain))
            qsub.writelines("#PBS -e /pandata/tlatrill/read/read_err_{0}_{1}\n".format(file_name, chain))
            qsub.writelines("#PBS -j oe\n")
            qsub.writelines("#PBS -W umask=022\n")
            qsub.writelines("#PBS -r n\n")
            qsub.writelines("#PBS -r n\n")
            qsub.writelines("TMP=/tmp/tlatrill$RANDOM\n")
            qsub.writelines("export TMPDIR=$TMP\n")
            for opt_name in ["mutsel", "mutselfreeomega"]:
                qsub.writelines("cd {}/pb_{}\n".format(data_path, opt_name))
                qsub.writelines("~/pbmpi2/data/readpb_mpi -x 100 -ss {0}\n".format(file_name))
                qsub.writelines("~/pbmpi2/data/readpb_mpi -x 100 -om {0}\n".format(file_name))
            qsub.writelines("rm -rf $TMP\n")
            qsub.writelines("rm {0}/qsub/{1}_{2}.pbs\n".format(data_path, file_name, chain))
            qsub.close()


print('Job completed')
