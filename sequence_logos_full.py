# Script title: sequence_logos
# Version title: full
# Short application: creates sequence logos of the full sequence while aligning them on there heptad positions of all df in 'E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta'
# Input files: csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py
# Output files: probability matrix, PWM and sequence logo plots of all the color_schemes given in color_schemes in 'E:\DENV_dimer_project\Output\sequence_logos_full_and_heptad'
# How to use: adjust color_scheme dictionary if necessary (look up other color_schemes of LOGO_MAKER) and run script. It automatically makes a sequence logo of all csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 23/09/2021

import pandas as pd
import itertools
import matplotlib.pyplot as plt
#plt.ion()
import logomaker as lm
import math
import os
from pathlib import Path

color_schemes = {"charge": "charge", "Chemistry": "chemistry", "hydrophobicity": "hydrophobicity",
                 "Skylign": "skylign_protein", "None": None} #change according to Logomaker package if you want other colorschemes

#finds first 7 heptad positions in register
def principal_period(s):
    i = (s+s).find(s, 1, -1)
    return None if i == -1 else s[:i]

#checks if CC is canonical (follows a strict abcdefg register)
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

#creates sequence logo of a certain colorscheme
def sequence_logo_full(CC_df, colorscheme = 'dimgray'):
    ind_max = 0
    for i in range(0, len(CC_df.poschain1)):
        if is_period(CC_df.poschain1[i]): #only includes CC if it is canonical
            letter = CC_df.poschain1[i][0]
            ind = "abcdefg".index(letter)
            ind = ind + len(CC_df.poschain1[i]) #checks width needed for the logo plot for one sequence
            if ind_max < ind: #maximizes the width needed of the plot to get the optimal plot size
                ind_max = ind
    xdim_matrix = ind_max

    # make a count matrix (row = heptad position and location in sequence, col = aminoacids)
    counttable = [["p", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]]
    for i in range(0, xdim_matrix): #builds empty matrix for logo
        listadd = [i, *itertools.repeat(0, 20)]
        counttable.append(listadd)

    col = ["pos", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    row = ["pos", *range(0, xdim_matrix)]

    #counts the amino acid occupation at the different positions of the sequence
    fig, ax = plt.subplots(1, 1, figsize=[4, 2])
    for chain in CC_df.chain1:
        hep_pos = CC_df[CC_df['chain1'] == chain]['poschain1']
        hep_pos = hep_pos.values[0]
        if is_period(hep_pos): #only includes a sequence if it is canonical
            hep_num = math.floor(len(hep_pos)/7) #checks from where the sequence should fill the logo
            sub = principal_period(hep_pos[:hep_num*7])
            if sub:
                n = "abcdefg".index(hep_pos[0])
                for i in range(0, len(chain)):
                    aa = chain[i]
                    pos = i + n
                    if aa != "X":
                        counttable[row.index(pos)][col.index(aa)] = counttable[row.index(pos)][col.index(aa)] + 1
            else:
                print(hep_pos[:hep_num*7], sub)

    #creates logo figure and probability matrix (PM) and position weighted matrix (PWM)
    df_logo = pd.DataFrame(counttable[1:], columns=counttable[0])
    del df_logo['p']
    df_logo_prob = lm.transform_matrix(df_logo, from_type='counts', to_type='probability') #transforms to PPM
    df_logo_weight = lm.transform_matrix(df_logo_prob, from_type='probability', to_type='weight')
    logo = lm.Logo(df_logo_prob, fade_probabilities=False, vpad=0.1, stack_order='small_on_top',
                   color_scheme=colorscheme, font_name='Arial', ax=ax) #boulds logo

    #creates list for ticks to ensure that the heptad position is printed on the x-axis
    xax = list(itertools.chain(*itertools.repeat(["a", "b", "c", "d", "e", "f", "g"], int(math.ceil((xdim_matrix)/7)))))
    erase = math.ceil((xdim_matrix) / 7)*7 - xdim_matrix #erases last heptad positions that are not filled
    if erase > 0:
        xax = xax[:-erase]
    logo.style_xticks(fmt='%d', anchor=0)
    logo.ax.set_xticklabels(x for x in xax) #sets ticks
    return logo, df_logo_prob, df_logo_weight


for filename in os.listdir('E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta'):
    print(filename[:-22])
    cdhit_df = pd.read_csv('E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta/' + filename, header=0)
    prob_matrix = sequence_logo_full(cdhit_df, colorscheme=None)[1] #gets probability matrix
    prob_matrix.to_csv('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/probability_matrix/seq_logo_full' + filename[:-22] + ".csv") #saves probability matrix
    plt.close()
    weight_matrix = sequence_logo_full(cdhit_df, colorscheme=None)[2] #gets weight matrix
    weight_matrix.to_csv('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/weight_matrix/seq_logo_full' + filename[:-22] + ".csv") #saves PWM
    plt.close()
    for color_scheme in color_schemes.keys(): #saves plot for all colors in color_schemes
        logo_out = sequence_logo_full(cdhit_df, colorscheme= color_schemes[color_scheme])
        logo = logo_out[0]
        logofig = plt.gcf()
        plt.close()
        plt.show()
        plt.close()
        Path('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/' + color_scheme).mkdir(parents=True, exist_ok=True)
        logofig.savefig('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/' + color_scheme + '/seq_logo_full_' + filename[:-22] + ".png") #saves plot logo