import os

path = "om_79_cds_mammals_no_pan_marsu_1"
tree = "mammals_no_homi_marsu.tree"
data = "/pandata/tlatrill/AdaptaPop/data"
# data = "./data"

os.system('rm ' + data + "/qsub/*")
for file in os.listdir(data + "/" + path):
    if file.endswith("filtered_NT.ali"):
        file_name = file[:-4]
        qsub = open(data + "/qsub/" + file_name + ".pbs", 'w')
        qsub.writelines("#!/bin/bash\n")
        qsub.writelines("#\n")
        qsub.writelines("#PBS -q q1day\n")
        qsub.writelines("#PBS -l nodes=1:ppn=8,mem=4gb\n")
        qsub.writelines("#PBS -o /pandata/tlatrill/out" + file_name + "\n")
        qsub.writelines("#PBS -e /pandata/tlatrill/err" + file_name + "\n")
        qsub.writelines("#PBS -j oe\n")
        qsub.writelines("#PBS -W umask=022\n")
        qsub.writelines("#PBS -r n\n")
        qsub.writelines("cp " + data + "/" + path + "/" + file + " ./" + file + "\n")
        qsub.writelines("cp " + data + "/" + tree + " ./" + tree + "\n")
        for chain in ["_1", "_2"]:
            command = "mpirun -n 8 ~/pbmpi/data/pb_mpi -f -x 1 200 -s -d ./" + file
            command += " -T ./" + tree + " -mutsel -dp ./" + file_name + chain + "\n"
            qsub.writelines(command)
            qsub.writelines("mv ./" + file_name + chain + ".* " + data + "/pb_output\n")
        qsub.close()

print('Job completed')
