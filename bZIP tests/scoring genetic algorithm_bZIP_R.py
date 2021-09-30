# Script title: scoring genetic algorithm_bZIPs
# Version title: Relative (R)
# Short application: scores bZIPs based on the scoring function relative scores
# Input files: csv file containing the position weighted matrix of the heptad lenght 24-28 PCCs, csv files of the interaction scoring matrices, register and sequence files of the bZIPs
# Output files: csv file containing the correlation of the scores and bZIPs
# How to use: change run name and run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 29/09/2021

import numpy as np
import pandas as pd
import os
import math

def listToString(s):
    str1 = ""
    return (str1.join(s))

#variables
run_name = "bZIP_Relative_Scores"

#input
Int_scoring_matrix = pd.read_csv("E:/DENV_dimer_project/Output/Interaction_volcanoplots/rel_prob_for_every_aa_int.csv", header=0, keep_default_na=False) #csv containing relative probability of interactions PCCvsAPCC
Pos_weight_matrix = pd.read_csv("E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/weight_matrix/seq_logo_heptadCCplus_PCC_length_24_28.csv", header=0, keep_default_na=False) #PWM for heptad positions in sequence logo

##scoring function
def get_up_down_helix(aa_code, heptad, orientation):
    pos_down = 0
    interaction_dict = {"I1": [], "I2_P": [], "I3_P": [], "I2_A" : [], "I3_A" : [], "I4": [], "I5": [], "I6": []}
    for aa in aa_code: #checks heptad position, possible interaction for every amino acid
        hp = heptad[pos_down]
        down_dict = {"a": ['I2_P', 'I2_A'], "b": "I0", "c": "I0", "d": ["I3_P", "I3_A"], "e": ["I1", "I5"], "f": "I0", "g": ["I4", "I6"]}
        int_code = down_dict[hp]
        if int_code != "I0":
            if orientation == "PCC":
                up_dict = {"I2_P": "a", "I3_A": "a", "I6": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d",
                           "I5": "d", "I4": "e", "I0": "f", "I1": "g"}  # parallel
                heptad_OR = heptad
                aa_code_OR = aa_code
            elif orientation == "APCC":
                up_dict = {"I2_P": "a", "I3_A": "a", "I5": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d",
                           "I6": "d", "I1": "e", "I0": "f", "I4": "g"}  # anti-parallel
                heptad_OR = heptad[::-1]
                aa_code_OR = aa_code[::-1]
            for int_c in int_code:
                up_hp = list(up_dict.values())[list(up_dict.keys()).index(int_c)]
                positions_up = [i for i, x in enumerate(heptad_OR) if x == up_hp]
                dis_small_4 = [abs(x - pos_down) < 4 for x in positions_up]
                if any(dis_small_4): #only forms an interaction if there is an interaction partner closer then 4 positions away.
                    index = dis_small_4.index(True)
                    pos_up = positions_up[index]
                    if int_c == "I2_P" or int_c == "I3_P" or int_c == "I1" or int_c == "I4":
                        interaction = ''.join(sorted(listToString(
                            (aa, aa_code_OR[pos_up]))))  # saves amino acids that interact symmetric interactions
                    else:
                        interaction = listToString(
                            (aa, aa_code_OR[pos_up]))  # saves amino acids that interact non-symmetric interactions
                    interaction_dict[int_c].append(interaction)  # puts interaction in dictionary
        pos_down = pos_down + 1
    return interaction_dict
i = 0
i_end = len(Pos_weight_matrix.index)
i_hep = 0
Pos_weight_matrix_list = []
while i < i_end:
    if i_hep == 7:
        i_hep = 0
    hep_pos = "abcdefg"[i_hep]
    pos_for_matrix = hep_pos + str(math.floor(i/7)+1)
    Pos_weight_matrix_list.append(pos_for_matrix)
    i += 1
    i_hep += 1
Pos_weight_matrix["pos"] = Pos_weight_matrix_list
def objective(aa_code, register):
    num_int_APCC = 0
    num_int_PCC = 0
    APCC1_score = 0
    PCC1_score = 0
    APCC2_score = 0
    PCC2_score = 0
    I1_APCC1_score = 0
    I2_P_APCC1_score = 0
    I2_A_APCC1_score = 0
    I3_P_APCC1_score = 0
    I3_A_APCC1_score = 0
    I4_APCC1_score = 0
    I5_APCC1_score = 0
    I6_APCC1_score = 0
    I1_PCC1_score = 0
    I2_P_PCC1_score = 0
    I2_A_PCC1_score = 0
    I3_P_PCC1_score = 0
    I3_A_PCC1_score = 0
    I4_PCC1_score = 0
    I5_PCC1_score = 0
    I6_PCC1_score = 0
    I1_APCC2_score = 0
    I2_P_APCC2_score = 0
    I2_A_APCC2_score = 0
    I3_P_APCC2_score = 0
    I3_A_APCC2_score = 0
    I4_APCC2_score = 0
    I5_APCC2_score = 0
    I6_APCC2_score = 0
    I1_PCC2_score = 0
    I2_P_PCC2_score = 0
    I2_A_PCC2_score = 0
    I3_P_PCC2_score = 0
    I3_A_PCC2_score = 0
    I4_PCC2_score = 0
    I5_PCC2_score = 0
    I6_PCC2_score = 0
    weight_scores = {"a": 0, "b": 0, "c": 0, "d" : 0, "e" : 0, "f": 0, "g": 0}
    for orientation_in in ["APCC", "PCC"]: #projects sequence on APCC and PCC structure
        int_list = get_up_down_helix(aa_code, register, orientation=orientation_in) #gets interactions based on orientation (APCC or PCC)
        if orientation_in == "APCC":
            cols = ["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_expected_asAPCC"] #relative probabilities that are taken into account
            for int in int_list.keys(): #looks at relative probabilities of an interaction observed in true APCCs
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                for AA_comb in int_list[int]:
                    num_int_APCC += 1
                    if int == "I1":
                        I1_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I1_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I2_P":
                        I2_P_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I2_P_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I3_P":
                        I3_P_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I3_P_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I2_A":
                        I2_A_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I2_A_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I3_A":
                        I3_A_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I3_A_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I4":
                        I4_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I4_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I5":
                        I5_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I5_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I6":
                        I6_APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I6_APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                    APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
        elif orientation_in == "PCC":
            cols = ["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_expected_asPCC"]
            for int in int_list.keys():
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                for AA_comb in int_list[int]:
                    num_int_PCC += 1
                    if int == "I1":
                        I1_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I1_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I2_P":
                        I2_P_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I2_P_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I3_P":
                        I3_P_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I3_P_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I2_A":
                        I2_A_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I2_A_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I3_A":
                        I3_A_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I3_A_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I4":
                        I4_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I4_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I5":
                        I5_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I5_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    elif int == "I6":
                        I6_PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                        I6_PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
                    PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])
                    PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])
    AA_total_weight = 0
    for i in range(0,len(aa_code)): #scores the aa sequence based on how well the heptad sequence logo is followed to avoid CC unlike amino acids to become more concentrated because of negative/positive APCC/PCC selection
        factor = 1
        pos_i = register[i] + str(1)
        AA_weight = Pos_weight_matrix[Pos_weight_matrix["pos"] == pos_i][aa_code[i]].values
        all_weights = Pos_weight_matrix.iloc[:, 1:][Pos_weight_matrix["pos"] == pos_i].values[0].tolist()
        all_weights = [abs(ele) for ele in all_weights]
        max_weight = max(all_weights)
        AA_weight_score = (AA_weight*factor)/max_weight #scales score to the max absolute weight possible
        weight_scores[pos_i[0]] += float(AA_weight_score)
        AA_total_weight += AA_weight_score
    AA_total_weight_score = AA_total_weight[0]/len(aa_code)
    I_scores = [I1_APCC1_score, I2_P_APCC1_score, I2_A_APCC1_score, I3_P_APCC1_score, I3_A_APCC1_score, I4_APCC1_score, I5_APCC1_score, I6_APCC1_score, I1_PCC1_score, I2_P_PCC1_score,
                    I2_A_PCC1_score, I3_P_PCC1_score, I3_A_PCC1_score, I4_PCC1_score, I5_PCC1_score, I6_PCC1_score, I1_APCC2_score, I2_P_APCC2_score, I2_A_APCC2_score, I3_P_APCC2_score, I3_A_APCC2_score, I4_APCC2_score, I5_APCC2_score, I6_APCC2_score, I1_PCC2_score, I2_P_PCC2_score,
                    I2_A_PCC2_score, I3_P_PCC2_score, I3_A_PCC2_score, I4_PCC2_score, I5_PCC2_score, I6_PCC2_score]
    APCC_score = (APCC1_score + APCC2_score)
    PCC_score = (PCC1_score + PCC2_score)
    score = 0
    return [APCC_score, PCC_score, score, AA_total_weight_score, I_scores, APCC1_score, PCC1_score, APCC2_score, PCC2_score, num_int_APCC, num_int_PCC]
def column(matrix, i):
    return [row[i] for row in matrix]


## running code + saving output as fasta file or
i_pop = 0
d = []
for filename in os.listdir('E:/DENV_dimer_project/Input/bZIP tests/sequence_files_half'):
    i_pop += 1
    #ATF6B = filename == "BACH2.txt"
    seq_file = open("E:/DENV_dimer_project/Input/bZIP tests/sequence_files_half/" + filename, "r")
    reg_file = open("E:/DENV_dimer_project/Input/bZIP tests/register_files_half/" + filename, "r")
    CC = seq_file.read().split("\n")[0]
    reg = reg_file.read().split("\n")[0]
    seq_file.close()
    reg_file.close()
    last_gen_score = objective(CC, reg)
    d.append(
        {
            'pdb' : filename[:-4],
            'aa_seq': CC,
            'APCC_score': last_gen_score[0],
            'PCC_score': last_gen_score[1],
            'Total_score': last_gen_score[2],
            'Weight_score': last_gen_score[3],
            'I1_APCC1_score' : last_gen_score[4][0],
            'I2_P_APCC1_score': last_gen_score[4][1],
            'I2_A_APCC1_score': last_gen_score[4][2],
            'I3_P_APCC1_score': last_gen_score[4][3],
            'I3_A_APCC1_score': last_gen_score[4][4],
            'I4_APCC1_score': last_gen_score[4][5],
            'I5_APCC1_score': last_gen_score[4][6],
            'I6_APCC1_score': last_gen_score[4][7],
            'I1_PCC1_score': last_gen_score[4][8],
            'I2_P_PCC1_score': last_gen_score[4][9],
            'I2_A_PCC1_score': last_gen_score[4][10],
            'I3_P_PCC1_score': last_gen_score[4][11],
            'I3_A_PCC1_score': last_gen_score[4][12],
            'I4_PCC1_score': last_gen_score[4][13],
            'I5_PCC1_score': last_gen_score[4][14],
            'I6_PCC1_score': last_gen_score[4][15],
            'I1_APCC2_score' : last_gen_score[4][16],
            'I2_P_APCC2_score': last_gen_score[4][17],
            'I2_A_APCC2_score': last_gen_score[4][18],
            'I3_P_APCC2_score': last_gen_score[4][19],
            'I3_A_APCC2_score': last_gen_score[4][20],
            'I4_APCC2_score': last_gen_score[4][21],
            'I5_APCC2_score': last_gen_score[4][22],
            'I6_APCC2_score': last_gen_score[4][23],
            'I1_PCC2_score': last_gen_score[4][24],
            'I2_P_PCC2_score': last_gen_score[4][25],
            'I2_A_PCC2_score': last_gen_score[4][26],
            'I3_P_PCC2_score': last_gen_score[4][27],
            'I3_A_PCC2_score': last_gen_score[4][28],
            'I4_PCC2_score': last_gen_score[4][29],
            'I5_PCC2_score': last_gen_score[4][30],
            'I6_PCC2_score': last_gen_score[4][31],
            'num_int_APCC': last_gen_score[9],
            'num_int_PCC': last_gen_score[10]
        }
    )

df = pd.DataFrame(d)
df.to_csv("E:/DENV_dimer_project/Output/bZIP tests/" + run_name + ".csv", index=False)

