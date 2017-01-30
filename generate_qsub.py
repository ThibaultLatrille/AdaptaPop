import os

path = "om_79_cds_mammals_no_pan_marsu"
tree = "mammals_no_homi_marsu.tree"

for file in os.listdir("./data/" + path):
    if file.endswith("filtered_NT.ali"):
        file_name = file[:-4]
        qsub = open("./qsub/" + file_name + ".pbs", 'w')
        qsub.writelines("#!/bin/bash\n")
        qsub.writelines("#\n")
        qsub.writelines("#PBS -q q1day\n")
        qsub.writelines("#PBS -l nodes=1:ppn=8,mem=4gb\n")
        qsub.writelines("#PBS -o out\n")
        qsub.writelines("#PBS -e err\n")
        qsub.writelines("#PBS -j oe\n")
        qsub.writelines("#PBS -W umask=022\n")
        qsub.writelines("#PBS -r n\n")
        qsub.writelines("ulimit ­v 4096000 ­T 8\n")
        qsub.writelines("cp /pandata/tlatrill/AdaptaPop/data/" + path + "/" + file + " ./" + file + "\n")
        qsub.writelines("cp /pandata/tlatrill/AdaptaPop/data/" + tree + " ./" + tree + "\n")
        command = "mpirun -n 8 /panhome/tlatrill/pbmpi/data/pb_mpi -f -x 1 10 -s -d ./" + file
        command += " -T ./" + tree + " -mutsel -dp ./" + file_name + "\n"
        qsub.writelines(command)
        qsub.writelines("cp ./" + file_name + ".* /pandata/tlatrill/AdaptaPop/data/pbmpi")
        qsub.close()

print('Job completed')
