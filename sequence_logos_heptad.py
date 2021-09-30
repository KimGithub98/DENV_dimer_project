# Script title: sequence_logos
# Version title: heptad
# Short application: creates sequence logos of the heptad positions of all df in 'E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta'
# Input files: csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py
# Output files: probability matrix, PWM and sequence logo plots of all the color_schemes given in color_schemes in 'E:\DENV_dimer_project\Output\sequence_logos_full_and_heptad'
# How to use: adjust color_scheme dictionary if necessary (look up other color_schemes of LOGO_MAKER) and run script. It automatically makes a sequence logo of all csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 23/09/2021

import pandas as pd
import itertools
import matplotlib.pyplot as plt
# plt.ion()
import logomaker as lm
import os
from pathlib import Path

color_schemes = {"charge": "charge", "Chemistry": "chemistry", "hydrophobicity": "hydrophobicity",
                 "Skylign": "skylign_protein", "None": None} #change according to Logomaker package if you want other colorschemes

#checks if a CC is canonical
def is_period(hep_pos):
    positions = ["a", "b", "c", "d", "e", "f", "g"]
    i= 0
    for pos in hep_pos:
        if i == 0:
            pos_num_theo = positions.index(pos)
            pos_num_true = positions.index(pos)
            period = (pos_num_theo == pos_num_true)
        else:
            pos_num_true = positions.index(pos)
            period = (pos_num_theo == pos_num_true)
        if pos == "g":
            pos_num_theo = 0
        else:
            pos_num_theo += 1
        i += 1
        if period is False:
            return False
    return True

#creates a sequence logo
def sequence_logo_heptad(CC_df, colorscheme='dimgray'):
    # make an empty count matrix (row = heptad position, col = aminoacids)
    counttable = [
        ["p", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
        ['a', *itertools.repeat(0, 20)],
        ['b', *itertools.repeat(0, 20)],
        ['c', *itertools.repeat(0, 20)],
        ['d', *itertools.repeat(0, 20)],
        ['e', *itertools.repeat(0, 20)],
        ['f', *itertools.repeat(0, 20)],
        ['g', *itertools.repeat(0, 20)]]

    col = ["pos", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
           "Y"]
    row = ["pos", "a", "b", "c", "d", "e", "f", "g"]

    #counts the amino acid occupation of all sequences
    fig, ax = plt.subplots(1, 1, figsize=[4, 2])
    for chain in CC_df.chain1:
        hep_pos = CC_df[CC_df['chain1'] == chain]['poschain1']
        hep_pos = hep_pos.values[0]
        if is_period(hep_pos): #excludes non-canonical helices
            for i in range(0, len(chain)):
                aa = chain[i]
                hp = hep_pos[i]
                if aa != "X": #excludes all unknown (X) amino acids
                    counttable[row.index(hp)][col.index(aa)] = counttable[row.index(hp)][col.index(aa)] + 1
        #else: (hep_pos)

    #creates probability matrix, position weighted matrix and logo plot
    df_logo = pd.DataFrame(counttable[1:], columns=counttable[0])
    df_logo['p'] = df_logo['p'].replace(["a", "b", "c", "d", "e", "f", "g"],
                                        [0, 1, 2, 3, 4, 5, 6])
    del df_logo['p']
    df_logo_prob = lm.transform_matrix(df_logo, from_type='counts', to_type='probability') #makes probability matrix
    df_logo_weight = lm.transform_matrix(df_logo_prob, from_type='probability', to_type='weight') #makes PMW

    logo = lm.Logo(df_logo_prob, fade_probabilities=False, vpad=0.1, stack_order='small_on_top',
                   color_scheme=colorscheme, font_name='Arial', ax=ax) #makes logo

    #adds heptad positions to the x-axis ticks
    logo.style_xticks(fmt='%d', anchor=0)
    logo.ax.set_xticklabels(x for x in ["a", "b", "c", "d", "e", "f", "g"])
    return logo, df_logo_prob, df_logo_weight


for filename in os.listdir('E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta'):
    cdhit_df = pd.read_csv('E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta/' + filename, header=0)
    prob_matrix = sequence_logo_heptad(cdhit_df, colorscheme=None)[1]  # gets probability matrix
    prob_matrix.to_csv('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/probability_matrix/seq_logo_heptad' + filename[:-22] + ".csv")  # saves probability matrix
    plt.close()
    weight_matrix = sequence_logo_heptad(cdhit_df, colorscheme=None)[2]  # gets weight matrix
    weight_matrix.to_csv('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/weight_matrix/seq_logo_heptad' + filename[:-22] + ".csv")  # saves PWM
    plt.close()
    for color_scheme in color_schemes.keys():
        logo = sequence_logo_heptad(cdhit_df, colorscheme= color_schemes[color_scheme])[0]
        logofig = plt.gcf()
        plt.close()
        plt.show()
        plt.close()
        Path('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/' + color_scheme).mkdir(parents=True, exist_ok=True)
        logofig.savefig('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/' + color_scheme + '/seq_logo_heptad_' + filename[:-22] + ".png")  # saves plot logo
