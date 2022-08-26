import argparse
import os
from pathlib import Path

def get_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', help='input raw_reads path can use gzip file', required=True)
    parser.add_argument('-d', help='Path of DB files folder', required=True)
    parser.add_argument('-o', help='output_folder_path', default="./mapping" ,required=False)

    args=parser.parse_args()
    return args

def main():
    global args
    args=get_args()


if __name__ == "__main__":
    main()

db_list = [str(x) for x in Path(args.d).iterdir()]
output_path = Path(args.o).resolve()
output_path.mkdir(parents=True, exist_ok=True)

for db in db_list:
    tmp_out = str(output_path / Path(db).stem) + ".bam"
    command = ["minimap2 -ax map-ont -N 100", db, args.i, "| samtools view -bS -F 0x4 -F 0x900 | samtools sort -o", tmp_out]
    os.system(' '.join(command))

# mappy로 긁거나 samtools를 돌려서 그냥바로 확인
# alignment coverage등을 계산해버리면 좋을 것 같음.
# 그래도 양심적으로 Coverage가 어느정도 되어야할테니.
