import pandas as pd
import os

if __name__ == '__main__':

    for folder in os.listdir():
        if not os.path.isdir(folder): continue
        for file in os.listdir(folder):
            path = os.path.join(folder, file)
            if not path.endswith("tsv.gz"): continue
            df = pd.read_csv(path, sep=",", low_memory=False)
            syn = len(df[(df["TYPE"] == "Syn") & (df["COUNT"] / df["SAMPLE_SIZE"] > 0.)])
            non_syn = len(df[(df["TYPE"] == "NonSyn") & (df["COUNT"] / df["SAMPLE_SIZE"] > 0.)])
            pct = 100 * syn / (syn + non_syn)
            print(f"File {path} contains {syn} synonymous mutations, accounting for {round(pct, 2)}%")
