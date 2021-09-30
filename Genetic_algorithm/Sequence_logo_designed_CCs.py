# Script title: Sequence logo designed CCs
# Version title: N/A
# Short application: creates sequence logo's of the final population of designed CC sequences
# Input files: fasta file of the final population
# Output files: .svg of the sequence logos
# How to use: Run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 30/09/2021

import pandas as pd
import itertools
import matplotlib.pyplot as plt
import logomaker as lm
import math
from Bio import SeqIO

def sequence_logo_full(fasta_file_sequences, colorscheme = 'dimgray'):
    xdim_matrix = 26
    counttable = [["p", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]]
    for i in range(0, xdim_matrix): #builds empty matrix for logo
        listadd = [i, *itertools.repeat(0, 20)]
        counttable.append(listadd)

    col = ["pos", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    row = ["pos", *range(0, xdim_matrix)]

    fig, ax = plt.subplots(1, 1, figsize=[4, 2])

    for sequence in fasta_file_sequences:
        for i in range(0, len(sequence.seq)):
                aa = sequence.seq[i]
                if i in [2,9,16,23]:
                    aa = "K"
                elif i in [5,12,19]:
                    aa = "Q"
                elif i in [6,13,20]:
                    aa = "A"

                counttable[row.index(i)][col.index(aa)] = counttable[row.index(i)][col.index(aa)] + 1

    print(counttable)
    df_logo = pd.DataFrame(counttable[1:], columns=counttable[0])
    del df_logo['p']
    df_logo_prob = lm.transform_matrix(df_logo, from_type='counts', to_type='probability') #transforms to PPM
    logo = lm.Logo(df_logo_prob, fade_probabilities=False, vpad=0.1, stack_order='small_on_top',
                   color_scheme=colorscheme, font_name='Arial', ax=ax) #boulds logo
    #creates list for ticks
    xax = list(itertools.chain(*itertools.repeat(["a", "b", "c", "d", "e", "f", "g"], int(math.ceil((xdim_matrix)/7)))))
    erase = math.ceil((xdim_matrix) / 7)*7 - xdim_matrix #erases last heptad positions that are not filled
    if erase > 0:
        xax = xax[:-erase]
    #xax.insert(0, "a")
    logo.style_xticks(fmt='%d', anchor=0)
    logo.ax.set_xticklabels(x for x in xax) #sets ticks
    return logo, df_logo_prob

# make a count matrix (row = heptad position, col = aminoacids)
color_schemes = {"charge": "charge", "Chemistry": "chemistry", "hydrophobicity": "hydrophobicity",
                 "Skylign": "skylign_protein", "None": None}

for color_scheme in color_schemes.keys(): #saves plot for all colors in color_schemes
    fasta_file_sequences = SeqIO.parse(open('E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/run_final.fasta'),'fasta')
    logo_out = sequence_logo_full(fasta_file_sequences, colorscheme= color_schemes[color_scheme])
    logo = logo_out[0]
    logofig = plt.gcf()
    plt.close()
    plt.show()
    logofig.savefig('E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/sequence_logos/' + color_scheme + ".png") #saves plot logo
    plt.close()