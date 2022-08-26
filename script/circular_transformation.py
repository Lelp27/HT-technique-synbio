"""
verified 22.02.15
Author : Seong Kun Bak (tjdrns27@kribb.re.kr)
"""

import argparse
import sys
import pysam
import mappy
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(
        prog="circular_transformation",
        description = "transform circular DNA reads to start in Same position.")
    parser.add_argument('-i', type=str, help="input_reads file, default value is stdin", default=sys.stdin)
    parser.add_argument('-r', type=str, help="reference fasta file to set as start position.", required=True)
    # -o flag말고 stdout 으로 빠지도록하자.
    parser.add_argument('-o', type=str, default=sys.stdout, help="output_dna_file_path")
    parser.add_argument(
        '-p', type=int, default=100, required=False,
        help="reference sequence's position in transformed sequence. default is 100 bp")
    parser.add_argument('-n', type=int, default=2, help="Iteration Number default i")
    parser.add_argument('--format', type=str, default='fastq', help='Sequence file foramt [fasta, fastq]')


    args = parser.parse_args()

    return args

def rotate(str, n):
    return str[n:] + str[:n]

def circulize(record, aligner, para, iter_num):
    # Recursive fucntion

    if iter_num <= 0:
        return (record)

    tmp = next(aligner.map(str(record.seq)))
    if tmp.is_primary:
        pass
    else:
        return ('continue')
    if tmp.strand == 1:
        tmp_dna = rotate(record, tmp.q_st - para)
    else:
        tmp_dna = rotate(record.reverse_complement(id=record), tmp.ctg_len - tmp.q_en - para)

    return (circulize(tmp_dna, aligner, para, iter_num-1))

def main():
    args=get_args()
    aligner = mappy.Aligner(args.r, preset="sr")

    if args.o != sys.stdout:
        output_dna = open(args.o, "w")
    else:
        output_dna = sys.stdout

    n = 0
    for record in SeqIO.parse(args.i, args.format):
        try:
            out_read = circulize(record, aligner, args.p, args.n).format('fastq')
            if out_read == "continue":
                continue
            output_dna.write(out_read)
        except:
            n+=1
            continue
    output_dna.close()

if __name__ == "__main__":
    main()
