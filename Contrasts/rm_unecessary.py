import os


def is_float(i):
    try:
        float(i)
        return True
    except:
        return False


if __name__ == '__main__':
    for species_f in os.listdir():
        if not os.path.isdir(species_f):
            continue
        for model_f in os.listdir(species_f):
            path_f = os.path.join(species_f, model_f)
            if not os.path.isdir(path_f):
                continue
            for file_name in sorted(os.listdir(path_f)):
                if ("ADAPTIVE" not in file_name) and ("NEARLY_NEUTRAL" not in file_name):
                    continue
                fpath = os.path.join(path_f, file_name)
                split = file_name.split("_")

                if int(split[0]) > 50:
                    os.remove(fpath)
                    print(f"{fpath} is removed")
                elif file_name.endswith(".messages") or file_name.endswith(".profile") or file_name.endswith(".err"):
                    os.remove(fpath)
                    print(f"{fpath} is removed")
                elif "ADAPTIVE" == split[1]:
                    i = split[2].replace(".", "_").split("_")[0]
                    if is_float(i):
                        os.remove(fpath)
                        print(f"{fpath} is removed")

