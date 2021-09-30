# Script title: bZIP preparation
# Version title: N/A
# Short application: Creates a .csv file with the Kd values of homodimer bZIPs from the Reinke et al. 2013 paper
# Input files: A sequence file (sequences, line = 15 and kd file, line = 14) downloaded from NCBI supplementary Supp_Data_1 and Supp_Data_2 respectivily
# Output files: csv file with interaction partners and Kd values & .txt file with fasta format (Uncommend line 53/54 to create a fasta file)
# How to use: Download the kd values/sequences files and run the script, no additional changes or input is needed. Uncommend line 53/54 to create a fasta file

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 27/09/21

import pandas as pd

kd_df = pd.read_csv("E:/DENV_dimer_project/Input/bZIP tests/KD_and_sequence_files/kd values.csv", header=0)
sequences = pd.read_csv("E:/DENV_dimer_project/Input/bZIP tests/KD_and_sequence_files/sequences.csv", header=0)

x = kd_df["Species1"] == "HS"
y = kd_df["Species2"] == "HS"
xy = [a and b for a, b in zip(x, y)]
kd_df = kd_df[xy]
list_CC = []
fasta = []

df_selec = [["Protein1", "Protein2", "Kd", "Length"]]

for i in range(0, len(kd_df)):
    CC = kd_df.iloc[i]
    protein1 = CC[0]
    protein2 = CC[1]
    kd = CC[2]
    if protein1 == protein2:
        seq1 = sequences["Protein sequence"][sequences["Name"] == protein1]
        seq2 = sequences["Protein sequence"][sequences["Name"] == protein2]
        if any(sequences["Name"] == protein1):
            if protein1 not in list_CC:
                fasta.append('>%s\n%s' % (protein1, seq1.values[0]))
                #print(CC)
                list_CC.append(protein1)
                df_selec.append([protein1, protein2, kd, len(seq1.values[0])])
        if any(sequences["Name"] == protein2):
            if protein2 not in list_CC:
                #print(protein2, seq2.values[0])
                fasta.append('>%s\n%s' % (protein2, seq2.values[0]))
                list_CC.append(protein2)
                df_selec.append([protein1, protein2, kd, len(seq2.values[0])])
        output_path = 'E:/DENV_dimer_project/Output/bZIP tests/KD_and_sequence_files/bZIPfasta.txt'

        #with open(output_path, 'w') as f:  #uncommend if you want to create a new fasta file, make sure there is no txt file with the same name that already exists, this will append this file
            #f.write('\n'.join(fasta))

df_selec = pd.DataFrame(df_selec[1:], columns= df_selec[0])
df_selec.to_csv("E:/DENV_dimer_project/Output/bZIP tests/KD_and_sequence_files/Kd.csv", index=False)
