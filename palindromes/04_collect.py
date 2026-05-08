import argparse
import os
from collections import defaultdict
import sys

data_dir = "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K"

samples = defaultdict(list)

def main():

    parser = argparse.ArgumentParser(prog="04_collect.py", description="collect results from palindrover")
    parser.add_argument('-p', "--path", help="path to directory with palindrover results", required=True)

    args = parser.parse_args()
    data_dir = args.path

    try:
        entries = os.listdir(data_dir)
        files = [entry for entry in entries if os.path.isfile(f"{data_dir}/{entry}")]
        files = [file for file in files if file.endswith(".pal")]
        for f in files:
            print(f"processing file {f}", file=sys.stderr)
            sample_name = f.split("_")[0]
            
            for line in open(f"{data_dir}/{f}", 'r'):
                line = line.strip() 

                if line.startswith("#"):
                    continue
            
                line = line.split("\t")
                print(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[9]}A\t{int(line[2]) - int(line[1])}\t{sample_name}")
                print(f"{line[0]}\t{line[5]}\t{line[6]}\t{line[9]}B\t{int(line[6]) - int(line[5])}\t{sample_name}")


    except Exception as e:
        print(f"Error while reading files in {data_dir}: {e}")

if __name__ == "__main__":
    main()