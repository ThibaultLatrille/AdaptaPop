def count_line(data_path, file_name):
    file = open('{0}/{1}'.format(data_path, file_name), 'r')
    counter = 0
    for line in file:
        if line[0] != "#":
            counter += 1
    file.close()
    return counter

data_path = "/mnt/sda1/AdaptaPop/data"

for file_name in ["Homo_sapiens_88_GRCh38_cds_all.fasta"]:
    print('{0}: {1} lines'.format(file_name, count_line(data_path, file_name)))
