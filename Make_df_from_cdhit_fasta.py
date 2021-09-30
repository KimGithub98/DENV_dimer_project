# Script title: Make_df_from_cdhit_fasta
# Version title: NA
# Short application: makes df from cdhit_output
# Input files: txt file with cdhit output, csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py
# Output files: cdhit_filtered df in Make_df_from_cdhit_fasta with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title>
# How to use: First run Unique_CC_fasta_for_cdhit, then remove redundant sequences as discribed in the protocol and save the fasta file this as text file in Input Make_df_from_cdhit_fasta. Then add the same lengthlist and filenames as used for Unique_CC_fasta_for_cdhit and run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 23/09/2021

import pandas as pd
import numpy as np
from Bio import SeqIO

length_list_in = (24,25,26,27,28)
files = ["CCplus_APCC", "CCplus_PCC"]


def make_cdhit_df(fasta_file_cdhit_output, CC_df, length_list, output_name):
    cdhit_filtered_df = pd.DataFrame(columns=["pdbID", "bundleID", "Molecule", "gene", "organism", "CC_length", "from_p", "to_p", "chain1", "chain2", "poschain1", "poschain2", "pub_Title"])
    i = 0
    Boullist = []
    for chain in CC_df['chain1']:
        lenght_boul = len(chain) in length_list
        Boullist.append(lenght_boul)
    length_df = CC_df[np.array(Boullist)]
    for sequence in fasta_file_cdhit_output:
        i = i + 1
        append_df = length_df[length_df['bundleID'] == sequence.id]
        append_df = append_df.values.tolist()[0]
        df_length = len(cdhit_filtered_df)
        cdhit_filtered_df.loc[df_length] = append_df
    cdhit_filtered_df.to_csv("E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta/" + output_name + ".csv", index=False)
    return cdhit_filtered_df

length_range = "_length_" + str(length_list_in[0]) + "_" + str(length_list_in[-1])
for filename in files:
    cdhit_filtered_output = SeqIO.parse(open('E:/DENV_dimer_project/Input/Make_df_from_cdhit_fasta/cdhit_fastas/' + filename + length_range + "_cdhit_out.txt"), 'fasta')
    CC_df_in = pd.read_csv("E:/DENV_dimer_project/Output/CCplus_database/Used_and_examples/" + filename + ".csv", header=0)
    cdhit_filtered_df = make_cdhit_df(fasta_file_cdhit_output=cdhit_filtered_output, CC_df=CC_df_in, length_list=length_list_in, output_name= filename + length_range + "_cdhit_filtered_df")

