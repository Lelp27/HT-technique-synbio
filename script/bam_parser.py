import argparse
import sys
import re
import pandas as pd
import pysam
import matplotlib.pyplot as plt

# Alignment_Start position = cigar_match.match(i.cigarstring).group()
# End_position = start_position + query_alignment_length
# No information about non-alignment sequence in bamfile
# so that with read_name and position information should find sequence.pysam

# query_cover = reference_length

def get_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', help='input bamfile path', default=sys.stdin, required=True)
    parser.add_argument('-o', help='png output path', default="./coverage.png" ,required=False)
    parser.add_argument('-t', type=str, help='Plot title', default=None ,required=False)
    parser.add_argument('-p', type=list, help='optional parameters : [bins, xlim_1, xlim_2]', default=None ,required=False)

    args=parser.parse_args()
    return args


def query_cover(query_df, reference):
    return(query_df["query_cover"]/reference[query_df["ref"]])

def bam_parser(input_bam):

    bamfile = pysam.AlignmentFile(input_bam)
    reference = {refs:lengths for lengths, refs in zip(bamfile.lengths, bamfile.references)}

    cigar_start_pos = re.compile("[0-9]+")
    query_list = []

    for record in bamfile:
        query_list.append([record.qname, [record.reference_name, record.is_reverse], record.flag,
                        int(cigar_start_pos.match(record.cigarstring).group()),
                        record.query_alignment_length, record.reference_length,
                        record.get_tag('NM')])
    df = pd.DataFrame(query_list)
    df.columns = ["qname", "ref", "flag", "start", "align_len", "query_cover", "NM"]
    df["coverage"] = df.apply(lambda x:query_cover(x, reference), axis=1)

    return df

def plot_coverage(bam_src, title=None, save_path=None, params = None):
    plt.clf()
    plt.hist(bam_src["coverage"], bins=200)
    plt.title(title)
    plt.savefig(save_path)

if __name__=="__main__":
    args = get_args()

    df = bam_parser(args.i)
    plot_coverage(df, title=args.t, save_path=args.o, params = args.p)