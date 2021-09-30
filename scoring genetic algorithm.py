import numpy as np
import pandas as pd
import os
import math

def listToString(s):
    str1 = ""
    return (str1.join(s))

#variables
AA_probability_matrix = pd.read_csv('e:/SAMcc/output/sequence logos/probability_matrix/seq_logo_heptadcdhit_CCplus_length24_28_APCC_HOMO_df.csv', header=0) #PPM of heptad logo
AA_probability_matrix.drop(['pos'],inplace=True, axis = 1)
print(AA_probability_matrix)
#print(AA_probability_matrix)
n_pop = 500 #size of population
register_pop = "defgabcdefgabcdefgabcdefga" #register of the population
length_CC_pop = 26 #number of aminoacids in the population
n_iter = 100 #number of generations/iterations
r_cross = 0.5 #probability of cross over
r_mut = 0.05 #probability of mutation
PCC_neg_sel_weight = 1.5 #weight of the negative selection score on PCC favored interactions
run_name = "run1"

## building initial population of AA based on the probabilities of amino acids in the heptad register
heptad = "abcdefg"
pop_per_pos = []
for i in range(0,length_CC_pop):
    hep_pos = register_pop[i]
    aa_prob = list(AA_probability_matrix.loc[heptad.index(hep_pos)].values)
    pop_pos = np.random.choice(list(AA_probability_matrix.columns.values), size= n_pop, replace = True, p = aa_prob) #selects aa depending on their probability on the heptad position
    pop_per_pos.append(pop_pos.tolist())

pop = []
for i in range(0, n_pop):
    list_pop = list(item[i] for item in pop_per_pos)
    STR = listToString(list_pop)
    pop.append(STR)

print(pop)

##scoring function
def get_up_down_helix4(aa_code, heptad, orientation):
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
                    interaction = listToString(sorted((aa, aa_code_OR[pos_up]))) #saves amino acids that interact
                    interaction_dict[int_c].append(interaction) #puts interaction in dictionary
        pos_down = pos_down + 1
    return interaction_dict

Int_scoring_matrix = pd.read_csv("E:/SAMcc/output/Interaction volcano plots_all_pairs/length24_28_/rel_prob_for_every_aa_int.csv", header=0) #csv containing relative probability of interactions PCCvsAPCC
#print(Int_scoring_matrix)
Int_NR_scoring_matrix = pd.read_csv("E:/SAMcc/output/Interaction volcano plots_all_pairs/Length24_28_/NR_prob_for_every_aa_int.csv", header=0) #csv containing non-relative probability of interactions PCCvsAPCC

Pos_weight_matrix = pd.read_csv("E:/SAMcc/output/sequence logos/weight_matrix/seq_logo_fullcdhit_CCplus_length24_28_APCC_HOMO_df.csv", header=0) #PWM for heptad positions in sequence logo
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

print(Pos_weight_matrix)

def objective(aa_code, register, PCC_neg_sel_weight):
    last_gen_scoreR = objectiveR(aa_code, register, PCC_neg_sel_weight)
    last_gen_scoreNR = objectiveNR(aa_code, register, PCC_neg_sel_weight)
    #- last_gen_scoreR[1][9] - last_gen_scoreR[1][11]
    #affinity_score = sum(last_gen_scoreR[1][16:23]) + sum(last_gen_scoreR[1][0:7]) - sum(last_gen_scoreR[1][8:15]) + last_gen_scoreR[0] * 10
    #extra_weight = 1
    #cc_score = -last_gen_scoreNR[1][24] * extra_weight - last_gen_scoreNR[1][29] * extra_weight - last_gen_scoreNR[1][8] * extra_weight - last_gen_scoreNR[1][13] * extra_weight + last_gen_scoreNR[1][2] * extra_weight + \
               #last_gen_scoreNR[1][2] * extra_weight + last_gen_scoreNR[1][4] * extra_weight - sum(last_gen_scoreNR[1][24:31]) - \
               #sum(last_gen_scoreR[1][24:31]) + sum(last_gen_scoreR[1][17:20]) - last_gen_scoreR[1][8] * extra_weight - \
               #last_gen_scoreR[1][13] * extra_weight - sum(last_gen_scoreR[1][14:15]) + last_gen_scoreR[1][2] * extra_weight + last_gen_scoreR[1][4] * extra_weight
    #total_score = affinity_score * cc_score
    #affinity_score = sum(last_gen_scoreR[1][16:21]) - last_gen_scoreR[1][9] - last_gen_scoreR[1][11] + sum(
    #    last_gen_scoreR[1][0:5]) - sum(last_gen_scoreR[1][8:13]) + last_gen_scoreR[0] * 10
    affinity_score = sum(last_gen_scoreR[1][0:5]) + sum(last_gen_scoreR[1][16:21]) - last_gen_scoreR[1][25] - last_gen_scoreR[1][27] + last_gen_scoreR[0]*3
    #specificity_score = sum(last_gen_scoreR[1][9:12]) - sum(last_gen_scoreR[1][8:15]) + sum(last_gen_scoreR[1][1:4])
    extra_weight = 2
    cc_score = last_gen_scoreNR[1][2] * extra_weight + last_gen_scoreNR[1][4] * extra_weight - \
               sum(last_gen_scoreR[1][8:15])/2 + sum(last_gen_scoreR[1][1:4]) + last_gen_scoreR[1][18] * extra_weight + last_gen_scoreR[1][17] + last_gen_scoreR[1][20] * extra_weight + \
               last_gen_scoreR[1][19]

    total_score = affinity_score + cc_score
    return [total_score, affinity_score, cc_score]
def objectiveR(aa_code, register, PCC_neg_sel_weight):
    num_int_APCC = 0
    num_int_PCC = 0
    APCC_score = 0
    PCC_score = 0
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
        int_list = get_up_down_helix4(aa_code, register, orientation=orientation_in) #gets interactions based on orientation (APCC or PCC)
        if orientation_in == "APCC":
            cols = ["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_shuffled_asAPCC"] #relative probabilities that are taken into account
            #, "CCplus_PCC_asAPCC" --> optional
            #["AA", "CCplus_APCC_asAPCC", "CCplus_PCC_asAPCC"]
            #["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_shuffled_asAPCC"]
            for int in int_list.keys(): #looks at relative probabilities of an interaction observed in true APCCs
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                int_weight = int_weight_APCC[int]
                #print(AA_scores)
                for AA_comb in int_list[int]:
                    num_int_APCC += 1
                    #print(np.append(AA_scores[AA_scores["AA"] == AA_comb][cols[1:3]].values[0], AA_scores[AA_scores["AA"] == AA_comb][cols[3]].values[0]*0.5))
                    #APCC_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:]].values[0]) # a high score means often observed in APCC compared to PCC
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
            cols = ["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_shuffled_asPCC"]
            #["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_shuffled_asPCC"]
            #["AA", "CCplus_APCC_asPCC", "CCplus_PCC_asPCC"]
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
        #if register[i] == "a" or register[i] == "d":
            #factor = 1.8 #forces following of the register logo more strictly (factor = n*times) on register position
        #elif register[i] == "g" or register[i] == "e":
            #factor = 0.5
        #else: factor = 0
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
        #print(AA_weight_score)
        weight_scores[pos_i[0]] += float(AA_weight_score)
        AA_total_weight += AA_weight_score
        if register[i] == "g":
            i_str += 1
    AA_total_weight_score = AA_total_weight[0]/len(aa_code)
    #AA_total_weight_score = 0
    # composition score
    total_comp_score = 0
    #for aa_comp in composition.columns:
    #    true_freq = aa_code.count(aa_comp)/18
    #    ideal_freq = composition.loc[0][composition.columns == aa_comp]
    #    comp_score = abs(ideal_freq[0]-true_freq)
    #    total_comp_score += comp_score
    I_scores = [I1_APCC1_score, I2_P_APCC1_score, I2_A_APCC1_score, I3_P_APCC1_score, I3_A_APCC1_score, I4_APCC1_score, I5_APCC1_score, I6_APCC1_score, I1_PCC1_score, I2_P_PCC1_score,
                    I2_A_PCC1_score, I3_P_PCC1_score, I3_A_PCC1_score, I4_PCC1_score, I5_PCC1_score, I6_PCC1_score, I1_APCC2_score, I2_P_APCC2_score, I2_A_APCC2_score, I3_P_APCC2_score, I3_A_APCC2_score, I4_APCC2_score, I5_APCC2_score, I6_APCC2_score, I1_PCC2_score, I2_P_PCC2_score,
                    I2_A_PCC2_score, I3_P_PCC2_score, I3_A_PCC2_score, I4_PCC2_score, I5_PCC2_score, I6_PCC2_score]
    #print(I_scores[20])
    #score = APCC_score/2 + PCC_score + AA_total_weight_score/2 - total_comp_score/5 #calculates total score, weights the CC total score (AA_total_weight_score) 2* more imporartant
    #print(APCC_score, AA_total_weight_score, total_comp_score, score)
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
        int_list = get_up_down_helix4(aa_code, register, orientation=orientation_in) #gets interactions based on orientation (APCC or PCC)
        if orientation_in == "APCC":
            cols = ["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_shuffled_asAPCC"] #relative probabilities that are taken into account
            #, "CCplus_PCC_asAPCC" --> optional
            #["AA", "CCplus_APCC_asAPCC", "CCplus_PCC_asAPCC"]
            #["AA", "CCplus_APCCvsPCC_asAPCC", "CCplus_shuffled_asAPCC"]
            for int in int_list.keys(): #looks at relative probabilities of an interaction observed in true APCCs
                AA_scores = Int_scoring_matrix[Int_scoring_matrix["Int"] == int][cols]
                int_weight = int_weight_APCC[int]
                #print(AA_scores)
                for AA_comb in int_list[int]:
                    #print(np.append(AA_scores[AA_scores["AA"] == AA_comb][cols[1:3]].values[0], AA_scores[AA_scores["AA"] == AA_comb][cols[3]].values[0]*0.5))
                    #APCC_score += sum(AA_scores[AA_scores["AA"] == AA_comb][cols[1:]].values[0]) # a high score means often observed in APCC compared to PCC
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
            cols = ["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_shuffled_asPCC"]
            #["AA", "CCplus_APCCvsPCC_asPCC", "CCplus_shuffled_asPCC"]
            #["AA", "CCplus_APCC_asPCC", "CCplus_PCC_asPCC"]
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
        #if register[i] == "a" or register[i] == "d":
            #factor = 1.8 #forces following of the register logo more strictly (factor = n*times) on register position
        #elif register[i] == "g" or register[i] == "e":
            #factor = 0.5
        #else: factor = 0
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
    #AA_total_weight_score = 0
    # composition score
    total_comp_score = 0
    #for aa_comp in composition.columns:
    #    true_freq = aa_code.count(aa_comp)/18
    #    ideal_freq = composition.loc[0][composition.columns == aa_comp]
    #    comp_score = abs(ideal_freq[0]-true_freq)
    #    total_comp_score += comp_score
    I_scores = [I1_APCC1_score, I2_P_APCC1_score, I2_A_APCC1_score, I3_P_APCC1_score, I3_A_APCC1_score, I4_APCC1_score, I5_APCC1_score, I6_APCC1_score, I1_PCC1_score, I2_P_PCC1_score,
                    I2_A_PCC1_score, I3_P_PCC1_score, I3_A_PCC1_score, I4_PCC1_score, I5_PCC1_score, I6_PCC1_score, I1_APCC2_score, I2_P_APCC2_score, I2_A_APCC2_score, I3_P_APCC2_score, I3_A_APCC2_score, I4_APCC2_score, I5_APCC2_score, I6_APCC2_score, I1_PCC2_score, I2_P_PCC2_score,
                    I2_A_PCC2_score, I3_P_PCC2_score, I3_A_PCC2_score, I4_PCC2_score, I5_PCC2_score, I6_PCC2_score]
    #print(I_scores[20])
    #print(APCC_score, AA_total_weight_score, total_comp_score, score)
    return [AA_total_weight_score, I_scores]



## selection
import random
# selects 4 candidates and does a 2x tournament to output 2 parents

## mutating data
# crossover two parents to create two children

## genetic algorithm
def column(matrix, i):
    return [row[i] for row in matrix]

## running code + saving output as fasta file or
from Bio import SeqIO

run_fasta = SeqIO.parse(open('e:/SAMcc/output/genetic_algorithm_output/' + run_name + ".fasta"), 'fasta')
print(run_fasta)


d = []
for CC in run_fasta:
    last_gen_score = objective(CC.seq, register_pop, PCC_neg_sel_weight)
    d.append(
        {
            'pdb' : CC.description,
            'aa_seq': CC.seq,
            'Total_score': last_gen_score[0],
            'affinity_score': last_gen_score[1],
            'specificity_score': last_gen_score[2],
            'logo_score': last_gen_score[3]
        }
    )

df = pd.DataFrame(d)
print(d)
df.to_csv("e:/SAMcc/output/genetic_algorithm_output/" + run_name + "_RE_final.csv", index=False)

