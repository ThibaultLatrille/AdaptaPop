import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sample', required=True, type=str, dest="sample", help="Call sample file")
    args = parser.parse_args()

    df = pd.read_csv(args.sample, sep="\t")
    df = df.groupby("pop").size().reset_index(name='counts')
    df.to_csv(args.sample + ".groupedby", index=False)
