# Script title: Unique_CC_fasta_for_cdhit
# Version title: NA
# Short application: creates fasta files of all unique sequences with a certain length specified in length list and saves in E:/DENV_dimer_project/Output/Unique_CC_fasta_for_cdhit/
# Input files: csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py
# Output files: fasta file of all unique sequences of a certain length
# How to use: Change length_list_in to the lengths (number of amino acids) that you want to extract. specify which files you want to use in files (file names in the E:/DENV_dimer_project/Output/CCplus_database/Used_and_examples/ directory)

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 23/09/2021

import pandas as pd
import numpy as np

##user input
length_list_in = (24,25,26,27,28)
files = ["CCplus_APCC", "CCplus_PCC"]


##saves fasta files
def save_unique_bundles_as_fasta(CC_df, length_list, output_name = 'PCC_length_fasta'):
    ## save unique bundles of selected length as fasta
    Boullist = []
    for chain in CC_df['chain1']:
        lenght_boul = len(chain) in length_list
        Boullist.append(lenght_boul)
    lenght_df = CC_df[np.array(Boullist)]
    print(lenght_df)
    fasta = []
    i = 0
    for sequence in np.unique(lenght_df['chain1']):
        i = i + 1
        id = lenght_df['bundleID'][lenght_df['chain1'] == sequence].values[0]
        fasta.append('>%s\n%s' % (id, sequence))

    output_path = 'E:/DENV_dimer_project/Output/Unique_CC_fasta_for_cdhit/' + output_name + ".fasta"
    print(output_path)
    with open(output_path, 'w') as f:
        f.write('\n'.join(fasta))

length_range = "_length_" + str(length_list_in[0]) + "_" + str(length_list_in[-1])
for filename in files:
    CC_df_in = pd.read_csv("E:/DENV_dimer_project/Output/CCplus_database/Used_and_examples/" + filename + ".csv", header=0)
    save_unique_bundles_as_fasta(CC_df= CC_df_in, length_list = length_list_in, output_name = filename + length_range + "_fasta_for_cdhit")
