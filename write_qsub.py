import os

ali_folder = "om_79_cds_mammals_no_pan_marsu_1"
tree = "mammals_no_homi_marsu.tree"
data_path = "/pandata/tlatrill/AdaptaPop/data"
# data_path = "./data"

ali_path = "{0}/{1}".format(data_path, ali_folder)
nbr_cpu = 4
os.system("rm {0}/qsub/*".format(data_path))
for file in os.listdir(ali_path):
    if file.endswith("filtered_NT.ali"):
        file_name = file[:-4]
        for option, opt_name in [("-globalomega", "globalomega"),
                                 ("-siteomega", "siteomega"),
                                 ("-mutsel -dp", "mutsel"),
                                 ("-mutsel -freeomega -dp", "mutselfreeomega")]:
            for chain in ["1", "2"]:
                qsub = open("{0}/qsub/{1}_{2}_{3}.pbs".format(data_path, file_name, opt_name, chain), 'w')
                qsub.writelines("#!/bin/bash\n")
                qsub.writelines("#\n")
                qsub.writelines("#PBS -q q1day\n")
                qsub.writelines("#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu))
                qsub.writelines("#PBS -o /pandata/tlatrill/out_err/out_{0}_{1}_{2}\n".format(file_name, opt_name, chain))
                qsub.writelines("#PBS -e /pandata/tlatrill/out_err/err_{1}_{1}_{2}\n".format(file_name, opt_name, chain))
                qsub.writelines("#PBS -j oe\n")
                qsub.writelines("#PBS -W umask=022\n")
                qsub.writelines("#PBS -r n\n")
                qsub.writelines("#PBS -r n\n")
                qsub.writelines("TMP=/tmp/tlatrill$RANDOM\n")
                qsub.writelines("export TMPDIR=$TMP\n")
                qsub.writelines("mkdir -p {0}/pb_{1}\n".format(data_path, opt_name))
                command = "mpirun -n {0} ~/pbmpi2/data/pb_mpi -f -s -x 1 400 {1}".format(nbr_cpu, option)
                command += " -d {0}/{1}".format(ali_path, file)
                command += " -T {0}/{1}".format(data_path, tree)
                command += " {0}/pb_{1}/{2}_{3}\n".format(data_path, opt_name, file_name, chain)
                qsub.writelines(command)
                qsub.writelines("rm -rf $TMP\n")
                qsub.writelines("rm {0}/qsub/{1}_{2}_{3}.pbs\n".format(data_path, file_name, opt_name, chain))
                qsub.close()


print('Job completed')
