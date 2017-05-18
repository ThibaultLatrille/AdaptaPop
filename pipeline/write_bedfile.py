from cds_libraries import build_dict_transcripts

# data = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

version = "88"
GRCh = "GRCh38"


def build_bedfile(_data_path, file_name, in_dict_transcripts):
    in_bedfile = open("{0}/{1}".format(_data_path,  file_name), 'w')
    cds_list = open('{0}/cds_homo_pan.out'.format(data_path), 'r')
    cds_list_header = cds_list.readline().replace('\n', '').split('\t')
    for line in cds_list:
        line_split = line.replace('\n', '').split('\t')
        cds_dict = dict(zip(cds_list_header, line_split))
        tr_id = cds_dict["HomoTr"]
        if in_dict_transcripts.get(tr_id):
            for line in in_dict_transcripts[tr_id].befile_lines():
                in_bedfile.write(line)

    in_bedfile.truncate()
    cds_list.close()
    in_bedfile.close()


dict_transcripts, _, _ = build_dict_transcripts(data_path, 'Homo_sapiens_{0}_{1}.gtf'.format(version, GRCh))
build_bedfile(data_path, '{0}_{1}_interval_cds.bed'.format(version, GRCh), dict_transcripts)

print('Job completed')
