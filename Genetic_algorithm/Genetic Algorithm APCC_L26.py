# Script title: Genetic Algorithm Final
# Version title: APCC_L26
# Short application: Final genetic Algorithm 
# Input files: csv file containing the AA probabilities of each heptad position (PM_heptadcdhit_CCplus_length24_28_APCC_HOMO_df.csv) csv file containing the position weighted matrix of the full lenght 24-28 APCCs (PMW_fullcdhit_CCplus_length24_28_APCC_HOMO_df.csv)
# Output files: csv file with the sequences in the final population and it's scores, fasta file with all the sequences
# How to use: Change variables if necessary and run genetic algorithm

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 29/09/2021

import numpy as np
import pandas as pd
import math

def listToString(s):
    str1 = ""
    return (str1.join(s))

##variables
n_pop = 500 #size of population
register_pop = "defgabcdefgabcdefgabcdefga" #register of the population
length_CC_pop = 26 #number of aminoacids in the population
n_iter = 40 #number of generations/iterations
r_cross = 0.5 #probability of cross over
r_mut = 0.05 #probability of mutation
PCC_neg_sel_weight = 1.5 #weight of the negative selection score on PCC favored interactions
run_name = "run1"

## building initial population of AA based on the probabilities of amino acids in the heptad register
AA_probability_matrix = pd.read_csv('E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/probability_matrix/seq_logo_heptadCCplus_APCC_length_24_28.csv', header=0) #PPM of heptad logo
AA_probability_matrix.drop(['pos'],inplace=True, axis = 1)

heptad = "abcdefg"
pop_per_pos = []
for i in range(0,length_CC_pop):
    hep_pos = register_pop[i]
    if hep_pos == "b" or hep_pos == "c":
        pop_pos = ['A' for i_pop in range(n_pop)]
    elif hep_pos == "f":
        pop_pos = ['Q' for i_pop in range(n_pop)]
    else:
        aa_prob = list(AA_probability_matrix.loc[heptad.index(hep_pos)].values)
        pop_pos = np.random.choice(list(AA_probability_matrix.columns.values), size=n_pop, replace=True, p=aa_prob)  # selects aa depending on their probability on the heptad position
        pop_pos = pop_pos.tolist()
    pop_per_pos.append(pop_pos)
pop = []
for i in range(0, n_pop):
    list_pop = list(item[i] for item in pop_per_pos)
    STR = listToString(list_pop)
    pop.append(STR)

print("starting population is: ")
print(pop)

##import Position weighted matrix and makes list of it
Pos_weight_matrix = pd.read_csv("E:/DENV_dimer_project/Output/sequence_logos_full_and_heptad/weight_matrix/seq_logo_fullCCplus_APCC_length_24_28.csv", header=0) #PWM for heptad positions in sequence logo
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

##scoring function
#imports relative and non-relative probability matrices
Int_scoring_matrix = pd.read_csv("E:/DENV_dimer_project/Output/Interaction_volcanoplots/rel_prob_for_every_aa_int.csv", header=0, keep_default_na=False) #csv containing relative probability of interactions PCCvsAPCC
Int_NR_scoring_matrix = pd.read_csv("E:/DENV_dimer_project/Output/Interaction_volcanoplots/NR_prob_for_every_aa_int.csv", header=0, keep_default_na=False) #csv containing non-relative probability of interactions PCCvsAPCC

#projection onto helical wheel and puts all interaction in a dictionary per interactiontype (I1-6)
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
                    #print(aa, aa_code_OR, pos_up)
                    if int_c == "I2_P" or int_c == "I3_P" or int_c == "I1" or int_c == "I4":
                        interaction = ''.join(sorted(listToString((aa, aa_code_OR[pos_up])))) #saves amino acids that interact symmetric interactions
                    else:
                        interaction = listToString((aa, aa_code_OR[pos_up])) #saves amino acids that interact non-symmetric interactions
                    interaction_dict[int_c].append(interaction) #puts interaction in dictionary
        pos_down = pos_down + 1
    return interaction_dict
def objective(aa_code, register, PCC_neg_sel_weight):
    last_gen_scoreR = objectiveR(aa_code, register, PCC_neg_sel_weight)
    last_gen_scoreNR = objectiveNR(aa_code, register, PCC_neg_sel_weight)
    affinity_score = sum(last_gen_scoreR[1][1:5]) + sum(last_gen_scoreR[1][17:21])
    specificity_score = last_gen_scoreR[1][0] + last_gen_scoreR[1][5] + last_gen_scoreR[1][16] + last_gen_scoreR[1][21] - last_gen_scoreR[1][25] - last_gen_scoreR[1][27] - sum(last_gen_scoreR[1][8:16])/2
    extra_weight = 2
    cc_score = last_gen_scoreNR[1][2] * extra_weight + last_gen_scoreNR[1][4] * extra_weight  + sum(last_gen_scoreR[1][1:5]) + last_gen_scoreR[1][18] * extra_weight + last_gen_scoreR[1][17] + last_gen_scoreR[1][20] * extra_weight + \
               last_gen_scoreR[1][19] + last_gen_scoreR[0]*3

    total_score = affinity_score + cc_score + specificity_score
    return [total_score, affinity_score, cc_score, specificity_score]
def objectiveR(aa_code, register, PCC_neg_sel_weight):
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
    int_weight_APCC = {"I1": 1, "I2_P": 1, "I3_P": 1, "I2_A": 1, "I3_A": 1, "I4": 1, "I5": 1, "I6": 1}
    int_weight_PCC = {"I1": 1, "I2_P": 1, "I3_P": 1, "I2_A": 1, "I3_A": 1, "I4": 1, "I5": 1, "I6": 1}
    for orientation_in in ["APCC", "PCC"]: #projects sequence on APCC and PCC structure
        int_list = get_up_down_helix(aa_code, register, orientation=orientation_in) #gets interactions based on orientation (APCC or PCC)
        if orientation_in == "APCC":
            cols = ["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_expected_asAPCC"] #relative probabilities that are taken into account
            for int in int_list.keys(): #looks at relative probabilities of an interaction observed in true APCCs
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                int_weight = int_weight_APCC[int]
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
                    APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])*int_weight
                    APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])*int_weight
        elif orientation_in == "PCC":
            cols = ["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_expected_asPCC"]
            for int in int_list.keys():
                int_weight = int_weight_PCC[int]
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
                    PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])*int_weight#/PCC_neg_sel_weight # a high score means often observed in PCC compared to APCC
                    PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])*int_weight#/PCC_neg_sel_weight
    AA_total_weight = 0
    i_str = 1
    for i in range(0,len(aa_code)): #scores the aa sequence based on how well the heptad sequence logo is followed to avoid CC unlike amino acids to become more concentrated because of negative/positive APCC/PCC selection
        factor = 1
        pos_i = register[i] + str(i_str)
        AA_weight = Pos_weight_matrix[Pos_weight_matrix["pos"] == pos_i][aa_code[i]].values
        all_weights = Pos_weight_matrix.iloc[:, 1:][Pos_weight_matrix["pos"] == pos_i].values[0].tolist()
        all_weights = [abs(ele) for ele in all_weights]
        max_weight = max(all_weights)
        if AA_weight >= 0:
            AA_weight_score = (AA_weight*factor)/(max_weight) #scales score to the max absolute weight possible
        else:
            AA_weight_score = -((abs(AA_weight) * factor)/(max_weight))  # scales score to the max absolute weight possible
        weight_scores[pos_i[0]] += float(AA_weight_score)
        AA_total_weight += AA_weight_score
        if register[i] == "g":
            i_str += 1
    AA_total_weight_score = AA_total_weight[0]/len(aa_code)
    I_scores = [I1_APCC1_score, I2_P_APCC1_score, I2_A_APCC1_score, I3_P_APCC1_score, I3_A_APCC1_score, I4_APCC1_score, I5_APCC1_score, I6_APCC1_score, I1_PCC1_score, I2_P_PCC1_score,
                    I2_A_PCC1_score, I3_P_PCC1_score, I3_A_PCC1_score, I4_PCC1_score, I5_PCC1_score, I6_PCC1_score, I1_APCC2_score, I2_P_APCC2_score, I2_A_APCC2_score, I3_P_APCC2_score, I3_A_APCC2_score, I4_APCC2_score, I5_APCC2_score, I6_APCC2_score, I1_PCC2_score, I2_P_PCC2_score,
                    I2_A_PCC2_score, I3_P_PCC2_score, I3_A_PCC2_score, I4_PCC2_score, I5_PCC2_score, I6_PCC2_score]
    return [AA_total_weight_score, I_scores, num_int_PCC, num_int_APCC]
def objectiveNR(aa_code, register, PCC_neg_sel_weight):
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
    int_weight_APCC = {"I1": 1, "I2_P": 1, "I3_P": 1, "I2_A": 1, "I3_A": 1, "I4": 1, "I5": 1, "I6": 1}
    int_weight_PCC = {"I1": 1, "I2_P": 1, "I3_P": 1, "I2_A": 1, "I3_A": 1, "I4": 1, "I5": 1, "I6": 1}
    for orientation_in in ["APCC", "PCC"]: #projects sequence on APCC and PCC structure
        int_list = get_up_down_helix(aa_code, register, orientation=orientation_in) #gets interactions based on orientation (APCC or PCC)
        if orientation_in == "APCC":
            cols = ["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_expected_asAPCC"] #relative probabilities that are taken into account
            for int in int_list.keys(): #looks at relative probabilities of an interaction observed in true APCCs
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                int_weight = int_weight_APCC[int]
                for AA_comb in int_list[int]:
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
                    APCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])*int_weight
                    APCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])*int_weight
        elif orientation_in == "PCC":
            cols = ["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_expected_asPCC"]
            for int in int_list.keys():
                int_weight = int_weight_PCC[int]
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                for AA_comb in int_list[int]:
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
                    PCC1_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:2]].values[0])*int_weight#/PCC_neg_sel_weight # a high score means often observed in PCC compared to APCC
                    PCC2_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[2:3]].values[0])*int_weight#/PCC_neg_sel_weight
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
    return [AA_total_weight_score, I_scores]

## perform selection
import random
# selects 4 candidates and does a 2x tournament to output 2 parents
def selection(pop, scores, k=2):
    selection_ix = random.sample(range(0, len(pop)), k = 4) #samples 4 aa sequences
    selection_ix1, selection_ix2 = selection_ix[0], selection_ix[1]
    if scores[selection_ix1] < scores[selection_ix2]: #tournament between seq 1 and 2
        selection_ix1 = selection_ix2
    selection_ix3, selection_ix4 = selection_ix[2], selection_ix[3]
    if scores[selection_ix3] < scores[selection_ix4]:  #tournament between seq 3 and 4
        selection_ix3 = selection_ix4
    return pop[selection_ix1], pop[selection_ix3]

## mutating data
# crossover two parents to create two children
def crossover(p1, p2, r_cross):
    # children are copies of parents by default
    c1, c2 = p1, p2
    # check for recombination
    if random.uniform(0,1) < r_cross:
        # select crossover point that is not on the end of the string
        pt = random.randint(1, len(p1)-2)
        # perform crossover
        c1 = p1[:pt] + p2[pt:]
        c2 = p2[:pt] + p1[pt:]
    return [c1, c2]

# mutation operator
def mutation(aa_code, r_mut):
    new_aa_code = ""
    for i in range(len(aa_code)):
        register = register_pop[i]
        if register in ["a", "d", "e", "g"]:
            random_number = random.uniform(0,1)
        elif register in ["b", "c", "f"]:
            random_number = 1
        aa_prob = list(AA_probability_matrix.loc[heptad.index(register)].values)
        if random_number < r_mut:
            new_aa = np.random.choice(list(AA_probability_matrix.columns.values), size=1, replace=True, p=aa_prob).tolist()[0] #mutates on heptad probability
            while new_aa == aa_code[i]:
                new_aa = np.random.choice(list(AA_probability_matrix.columns.values), size=1, replace=True,p=aa_prob).tolist()[0]
            new_aa_code = new_aa_code + new_aa
        else: new_aa_code = new_aa_code + aa_code[i]
    return new_aa_code

## genetic algorithm
def column(matrix, i):
    return [row[i] for row in matrix]

def genetic_algorithm(pop, register_pop, n_pop, n_iter, r_cross, r_mut, PCC_neg_sel_weight):
    run_data = [["gen", "A_total_score", "A_affinity_score", "A_cc_score", "A_specificity_score", "best_total_score", "best_affinity_score", "best_cc_score", "best_specificity_score"]]
    best, best_score, best_affinity, best_specificity, best_logo = 0, 0, 0, 0, 0
    for gen in range(n_iter):
        all_scores = [objective(aa_code_in, register_pop, PCC_neg_sel_weight) for aa_code_in in pop]
        scores_affinity = column(all_scores, 1)
        scores_cc = column(all_scores, 2)
        scores = column(all_scores, 0)
        scores_specificity = column(all_scores, 3)
        for i in range(0, n_pop):
            if scores[i] > best_score:
                best, best_score, best_affinity, best_cc, best_specificity = pop[i], scores[i], scores_affinity[i], scores_cc[i], scores_specificity[i]
        print(">%d, new best %s = %.3f" % (gen,  best, best_score))
        print("best affinity score: ", best_affinity, ", best cc score: ", best_cc, ", best affinity score: ", best_specificity)
        run_data.append([gen, sum(scores)/len(scores), sum(scores_affinity)/len(scores_affinity), sum(scores_cc)/len(scores_cc), sum(scores_specificity)/len(scores_specificity), best_score, best_affinity, best_cc, best_specificity])
        selected = [selection(pop, scores) for _ in range(int(n_pop/2))]
        children = list()
        for i in range(0, int(n_pop/2)):
            p1, p2 = selected[i]
            for c in crossover(p1, p2, r_cross):
                child = mutation(c, r_mut)
                children.append(child)
        pop = children
    return pop, run_data

## running code + saving output as fasta file or
def save_unique_CCs_as_fasta(seq_list, filename):
    fasta = []
    i = 0
    for sequence in np.unique(seq_list):
        i = i + 1
        id = "CC" + str(i)
        fasta.append('>%s\n%s' % (id, sequence))
    output_path = 'e:/SAMcc/output/genetic_algorithm_output/' + filename + ".fasta"
    with open(output_path, 'w') as f:
        f.write('\n'.join(fasta))

out = genetic_algorithm(pop = pop, register_pop = register_pop, n_pop = n_pop, n_iter = n_iter, r_cross= r_cross, r_mut= r_mut, PCC_neg_sel_weight = PCC_neg_sel_weight)
print(out)
last_gen = out[0]
last_gen_score = [objective(aa_code_in, register_pop, PCC_neg_sel_weight) for aa_code_in in last_gen]
print(last_gen_score[0])
d = []
for i in range(0, len(last_gen)):
    d.append(
        {
            'aa_seq': last_gen[i],
            'total_score': last_gen_score[i][0],
            'affinity_score': last_gen_score[i][1],
            'cc_score': last_gen_score[i][2]
        }
    )

df = pd.DataFrame(d)
print(d)
df.to_csv("E:/DENV_dimer_project/Output/Genetic Algorithm/" + run_name, index=False)
save_unique_CCs_as_fasta(last_gen, run_name)
["gen", "A_total_score", "A_affinity_score", "A_cc_score", "A_specificity_score", "best_total_score", "best_affinity_score", "best_cc_score", "best_specificity_score"]

##visualize run data
import matplotlib.pyplot as plt
run_data = out[1]
run_data = pd.DataFrame(run_data[1:], columns= run_data[0])
plt.plot(run_data.gen, run_data.A_total_score, label = "Mean score")
plt.plot(run_data.gen, run_data.A_affinity_score, label = "Mean affinity score")
plt.plot(run_data.gen, run_data.A_cc_score, label = "Mean cc score")
plt.plot(run_data.gen, run_data.A_specificity_score, label = "Mean specificity score")
plt.plot(run_data.gen, run_data.best_total_score, label = "best total score")
plt.plot(run_data.gen, run_data.best_affinity_score, label = "best affinity score")
plt.plot(run_data.gen, run_data.best_cc_score, label = "best CC score")
plt.plot(run_data.gen, run_data.best_specificity_score, label = "best specificity score")
plt.xlabel('generation')
plt.legend()
plt.show()



