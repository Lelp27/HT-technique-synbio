"""
DNA Part Library seq DB 제작하기

# verified 22.02.17
Author : Seong Kun Bak

1. 합칠 DNA 조각의 수 선택
2. 각 조각 별 Library의 수 선택 (각 조각이 Library form 인지 or not)
3. Junction on or not (Junction 서열 선택)
4. Assemble and export to fasta
"""

# argument 없이 우선 구현하기.
# input을 fasta를 받을 것인가 str sequence or xlsx 를 받을 것인가.
import pandas as pd
from Bio import SeqIO

def df_to_fasta(name, seq, save_path):
    f = open(save_path, 'w')

    for i, j in zip(name, seq):
        fasta = '>' + i.strip('\n') + '\n' + j.strip('\n') + '\n'
        f.writelines(fasta)
    f.close()

def rev_com(seq):  #reverse complement for string.
    seq = seq[::-1]
    table = str.maketrans('ATGC','TACG')
    return seq.translate(table)

def na_to_blank(df):
    if type(df) == list:
        for i in df:
            i.fillna('', inplace=True)
    else:
        df.fillna('', inplace=True)

db_path = "/mnt/kun/Part_DB.xlsx"
pro_df = pd.read_excel(db_path, sheet_name = "pro", usecols=range(0,3))
rbs_df = pd.read_excel(db_path, sheet_name = "rbs", usecols=range(0,3))
ter_df = pd.read_excel(db_path, sheet_name = "ter", usecols=range(0,3))
cds_df = pd.read_excel(db_path, sheet_name = "cds", usecols=range(0,3))
linker_df = pd.read_excel(db_path, sheet_name="overhang", usecols=range(0,3))

na_to_blank([pro_df, rbs_df, ter_df, cds_df, linker_df])

#tmp = 'p21'

#pro_lib = [tmp, 'p22']
pro_lib = ['p9','p12','p21']
rbs_lib = ['r11','r22']
ter_lib = ['t25']
cds_lib = ["sfGFP"]
linker_lib = ["O2", "O3", "O4"]
#cds_lib = ["tphA1", "tphA2", "tphA3", "tphB"]

#save_path = '/mnt/kun/raw_data/DB/TCC/220325/p117_c32_tphA3/' + tmp + '.fasta'
save_path = '/mnt/kun/gfp6/gfp6.fasta'

pro_df = pro_df[pro_df["No"].isin(pro_lib)]
rbs_df = rbs_df[rbs_df["No"].isin(rbs_lib)]
ter_df = ter_df[ter_df["No"].isin(ter_lib)]
cds_df = cds_df[cds_df["Name"].isin(cds_lib)]
# linker_lib의 order가 그대로 유지되지 않음.
linker_df = linker_df[linker_df["No"].isin(linker_lib)]
linker_df['No'] = pd.Categorical(linker_df['No'], linker_lib)
linker_df.sort_values('No', inplace=True)
linker_df = linker_df["Sequence"].values

#empty Dataframe
final_df = pd.DataFrame({"Name" : [], "Seq" : []})
final_list = []
for i1 in pro_df.iterrows():
    for i2 in rbs_df.iterrows():
        for i3 in cds_df.iterrows():
            for i4 in ter_df.iterrows():
                name = "_".join([
                    i1[1]["No"], i2[1]["No"], i3[1]["Name"], i4[1]["No"]])
                seq = "aatcgtctcaaaaggcct" + i1[1]["Sequence"] + linker_df[0] + i2[1]["Sequence"] + linker_df[1] + i3 [1]["Sequence"] + linker_df[2] + i4[1]["Sequence"]
                seq = seq.upper()
                final_list.append(pd.DataFrame({"Name":[name], "Seq":[seq]}))
                final_df = pd.concat(final_list, ignore_index=True)
                #final_df = final_df.append({"Name":name, "Seq":seq}, ignore_index=True)

df_to_fasta(name = final_df["Name"].values, seq = final_df["Seq"].values, save_path=save_path)