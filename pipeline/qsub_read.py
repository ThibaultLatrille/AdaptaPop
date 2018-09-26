import os
from subprocess import run

ali_folder = "pb_cleaned_mutsel"
data_path = "/mnt/sda1/AdaptaPop/data"
pb_mpi_path = "/home/thibault/Tools/pbmpi2/data"
# data_path = "./data"

ali_path = "{0}/{1}".format(data_path, ali_folder)
for file in os.listdir(ali_path):
    if file.endswith(".chain"):
        file_name = file[:-len(".chain")]
        qsub_path = "{0}/{1}.pbs".format(data_path, file_name)
        qsub = open(qsub_path, 'w')
        qsub.writelines("#!/bin/bash\n")
        qsub.writelines("#\n")
        qsub.writelines("#PBS -q q1day\n")
        qsub.writelines("#PBS -l nodes=1:ppn=1,mem=1gb\n")
        qsub.writelines("#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(file_name))
        qsub.writelines("#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(file_name))
        qsub.writelines("#PBS -j oe\n")
        qsub.writelines("#PBS -W umask=022\n")
        qsub.writelines("#PBS -r n\n")
        qsub.writelines("#PBS -r n\n")
        qsub.writelines("TMP=/tmp/tlatrill$RANDOM\n")
        qsub.writelines("export TMPDIR=$TMP\n")
        for opt_name in ["mutsel"]:
            qsub.writelines("cd {0}\n".format(ali_path))
            qsub.writelines("mpirun -n 8 {0}/readpb_mpi -x 100 -ss {1}\n".format(pb_mpi_path, file_name))
            # qsub.writelines("{0}/readpb_mpi -x 100 -om {1}\n".format(pb_mpi_path, file_name))
        qsub.writelines("rm -rf $TMP\n")
        qsub.writelines("rm {0}\n".format(qsub_path))
        qsub.close()

        run("sh {0}".format(qsub_path), shell=True)


print('Job completed')
