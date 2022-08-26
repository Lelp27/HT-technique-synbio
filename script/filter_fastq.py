import argparse
import sys
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description = "Filter fastq reads with query_name files")
    parser.add_argument('-i', type=str, default=sys.stdin, help="input DNA file_path want to filter")
    parser.add_argument('-q', type=str, help="input query_name file path", required=True)
    # not available
    parser.add_argument('-o', default=sys.stdout, help="output file Path")
    
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    input_dna = SeqIO.parse(args.i, "fastq")
    with open(args.q, "rt") as f:
        qname = f.readlines()
    qname = [i.strip() for i in qname]

    if args.o != sys.stdout:
        output_dna = open(args.o, "w")
    else:
        output_dna = sys.stdout

    for record in input_dna:
        if record.name in qname:
            output_dna.write(record.format("fastq"))
        else:
            pass
    output_dna.close()

if __name__ == "__main__":
    main()