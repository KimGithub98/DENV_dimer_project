import pandas as pd
import os

#print(os.getcwd())
cdhit_df = pd.read_csv('e:/SAMcc/output/cdhit_dfs/cdhit_CCplus_length14_18_APCC_HOMO_df_selection.csv', header=0)

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'UNK': 'X'}
inv_d = {v: k for k, v in d.items()}

for CC in cdhit_df["bundleID"]:
    list_file = open("mutation_list.txt", "a")
    length = cdhit_df[cdhit_df["bundleID"] == CC]["CC_length"]
    length = length.values[0]
    if length == 17:
        print(CC[0:4])
        i = 0
        aa_seq = cdhit_df[cdhit_df["bundleID"] == CC]["chain1"]
        aa_seq = aa_seq.values[0]
        #adds to mutate list original pdb , chain, resnumber, from AA, to AA
        for aa in aa_seq:
            i += 1
            for chain in ("A", "B"):
                if chain == "B":
                    letter = i + 18
                if chain == "A":
                    letter = i
                new_line = CC[0:4] + ".pdb" + " " + chain + " " + str(letter) + " " + "GLY" + " " + inv_d[aa] + " " + "\n"
                #print(new_line)
                list_file.write(new_line)
    list_file.close()