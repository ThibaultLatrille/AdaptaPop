#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    ali_file = open(args.input, 'r')
    ali_file.readline()
    outfile = open(args.output, "w")
    outfile.write("\n".join(
        [">{0}\n{1}".format(*line.replace("  ", " ").replace("\t", " ").strip().split(" ")) for line in
         ali_file.readlines()]))
    outfile.close()
