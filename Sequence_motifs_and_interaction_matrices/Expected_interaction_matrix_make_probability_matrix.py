# Script title: Expected_interaction_matrix
# Version title: Make_probability_matrix
# Short application: Makes a probability matrix of expected interaction values based on the amino acid frequencies and length and starting register of the CCs
# Input files: csv file with csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title> created by CCplus prepare.py, csv file describing the interaction_bias for each length and starting position
# Output files: csv files with the expected interactions
# How to use: Change interaction bias, CCplus__cdhit_filtered_df and frequencies of amino acids if neccasary and run script. It automatically creates a file for PCC and APCC

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 24/09/2021

import pandas as pd
import itertools

#input
interaction_bias = pd.read_csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/interaction_bias_24_28.csv", header=0)
CCplus = pd.read_csv('E:/DENV_dimer_project/Output/Make_df_from_cdhit_fasta/CCplus_APCC_length_24_28_cdhit_filtered_df.csv', header=0)

#checks if a register is in period
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

#function to get all possible combinations with replacement (a.i. AC and CA)
def permutations_with_replacement(m, n):
    for i in itertools.product(m, repeat=n):
        yield i

#creates empty interaction count matrix (rows = amino acid combination, cols = Interaction types)
def empty_interaction_count_matrix():
    ## make empty interaction count matrix
    aalist = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    interactions = permutations_with_replacement(aalist, 2)
    interactions = [''.join(i) for i in interactions]
    interactions.insert(0, "pos")
    int_count_matrix = [interactions]
    int_list = ["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"]
    for int in int_list:
        listadd = [int, *itertools.repeat(0, 400)]
        int_count_matrix.append(listadd)
    return int_count_matrix

#creates lists with bias for each interaction for a CC depending on length and starging position
I4 = []
I1 = []
I2_P_pcc = []
I3_P_pcc = []
I2_P_apcc = []
I3_P_apcc = []
for i in range(0,len(CCplus.poschain1)):
    if is_period(CCplus.poschain1[i]):
        bias_code = str(len(CCplus.poschain1[i])) + str(CCplus.poschain1[i][0]) #<length><starting heptad position> code to search the interaction_bias matrix
        bias_for_CC = interaction_bias[interaction_bias.Coiled_coil == bias_code] #looks for percentage of bias for that length and start register
        # appends the bias percentage to a list for every interaction
        I4.append(bias_for_CC.I4.values[0])
        I1.append(bias_for_CC.I1.values[0])
        I3_P_apcc.append(bias_for_CC.I3_P_apcc.values[0])
        I2_P_apcc.append(bias_for_CC.I2_P_apcc.values[0])
        I3_P_pcc.append(bias_for_CC.I3_P_pcc.values[0])
        I2_P_pcc.append(bias_for_CC.I2_P_pcc.values[0])

#calculates average bias in the datasets given at input
CCplus_bias_I1 = sum(I1)/len(I1)
CCplus_bias_I4 = sum(I4)/len(I4)
CCplus_bias_I2_P_apcc = sum(I2_P_apcc)/len(I2_P_apcc)
CCplus_bias_I3_P_apcc = sum(I3_P_apcc)/len(I3_P_apcc)
CCplus_bias_I2_P_pcc = sum(I2_P_pcc)/len(I2_P_pcc)
CCplus_bias_I3_P_pcc = sum(I3_P_pcc)/len(I3_P_pcc)
CCplus_bias_allI = [CCplus_bias_I1, CCplus_bias_I2_P_pcc, CCplus_bias_I3_P_pcc, CCplus_bias_I2_P_apcc, CCplus_bias_I3_P_apcc, CCplus_bias_I4]


#calculates expected interactions depending on the frequencies in CCs
frequencies = {"A": 0.0847, "C": 0.00745, "D": 0.0456, "E": 0.139, "F": 0.0289, "G": 0.0158, "H": 0.0233, "I": 0.0661, "K": 0.0745, "L": 0.161, "M": 0.027, "N": 0.0307, "P": 0.00372, "Q": 0.0773, "R": 0.0717, "S": 0.0512, "T": 0.0289, "V": 0.0447, "W": 0.00372, "Y": 0.0149}
## make empty interaction count matrix
interactions = permutations_with_replacement(frequencies.keys(), 2)
interactions = [''.join(i) for i in interactions]

## make empty interaction count matrix
sumP = 0
down_dict = {"I2_P": "a", "I3_A": "d", "I5": "e", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "a", "I6": "g", "I1": "e", "I0": "f", "I4": "g"}  # anti-parallel
hep_res_to_pos = {"a": 0, "b": 1, "c":2, "d":3, "e":4, "f":5, "g":6}
int_codes = ["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"]
frequency_matrix_APCC = empty_interaction_count_matrix()
frequency_matrix_PCC = empty_interaction_count_matrix()

#calculates probability interaction matrix based on the frequency of the amino acids and bias
for orientation in ["APCC", "PCC"]:
    if orientation == "APCC":
        up_dict = {"I2_P": "a", "I3_A": "a", "I5": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I6": "d", "I1": "e", "I0": "f", "I4": "g"}  # anti-parallel
        for int_code in int_codes:
            hep_res1 = hep_res_to_pos[down_dict[int_code]]
            hep_res2 = hep_res_to_pos[up_dict[int_code]]
            frequencies_res1 = frequencies
            frequencies_res2 = frequencies
            row_index = int_codes.index(int_code) + 1
            if int_code == "I4":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1]*frequencies_res2[res2]) * (1-CCplus_bias_I4)
                        intP_bias = (frequencies_res1[res1] * frequencies_res2[res2])
                    elif res1 == res2:
                        intP = frequencies_res1[res1] * frequencies_res2[res2] * (1-CCplus_bias_I4) + frequencies_res1[res1] * (CCplus_bias_I4)
                        intP_bias = frequencies_res1[res1] * frequencies_res2[res2]
                    col_num = frequency_matrix_APCC[0].index(''.join(sorted(int)))
                    frequency_matrix_APCC[row_index][col_num] += intP
            elif int_code == "I1":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1]*frequencies_res2[res2]) * (1-CCplus_bias_I1)
                    elif res1 == res2:
                        #print(res1, frequencies[res1])
                        intP = (frequencies_res1[res1] * frequencies_res2[res2] * (1-CCplus_bias_I1)) + (frequencies_res1[res1] * (CCplus_bias_I1))
                        intP_bias = frequencies_res1[res1] * frequencies_res2[res2]
                        #print(intP, intP_bias)
                    col_num = frequency_matrix_APCC[0].index(''.join(sorted(int)))
                    frequency_matrix_APCC[row_index][col_num] += intP
            elif int_code == "I2_P":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1]*frequencies_res2[res2]) * (1-CCplus_bias_I2_P_apcc)
                    elif res1 == res2:
                        #print(res1, frequencies[res1])
                        intP = (frequencies_res1[res1] * frequencies_res2[res2] * (1-CCplus_bias_I2_P_apcc)) + (frequencies_res1[res1] * (CCplus_bias_I2_P_apcc))
                        intP_bias = frequencies_res1[res1] * frequencies_res2[res2]
                        #print(intP, intP_bias)
                    col_num = frequency_matrix_APCC[0].index(''.join(sorted(int)))
                    frequency_matrix_APCC[row_index][col_num] += intP
            elif int_code == "I3_P":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1] * frequencies_res2[res2]) * (
                                    1 - CCplus_bias_I3_P_apcc)
                    elif res1 == res2:
                        # print(res1, frequencies[res1])
                        intP = (frequencies_res1[res1] * frequencies_res2[res2] * (1 - CCplus_bias_I3_P_apcc)) + (
                                    frequencies_res1[res1] * (CCplus_bias_I3_P_apcc))
                        intP_bias = frequencies_res1[res1] * frequencies[res2]
                        # print(intP, intP_bias)
                    col_num = frequency_matrix_APCC[0].index(''.join(sorted(int)))
                    frequency_matrix_APCC[row_index][col_num] += intP
            else:
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = frequencies_res1[res1]*frequencies_res2[res2]
                    elif res1 == res2:
                        intP = frequencies_res1[res1] * frequencies_res2[res2]
                    col_num = frequency_matrix_APCC[0].index(int)
                    frequency_matrix_APCC[row_index][col_num] += intP
    if orientation == "PCC":
        up_dict = {"I2_P": "a", "I3_A": "a", "I6": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I5": "d", "I4": "e", "I0": "f", "I1": "g"}  # parallel
        for int_code in int_codes:
            hep_res1 = hep_res_to_pos[down_dict[int_code]]
            hep_res2 = hep_res_to_pos[up_dict[int_code]]
            frequencies_res1 = frequencies
            frequencies_res2 = frequencies
            row_index = int_codes.index(int_code) + 1
            if int_code == "I2_P":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1]*frequencies_res2[res2]) * (1-CCplus_bias_I2_P_pcc)
                        intP_bias = (frequencies_res1[res1] * frequencies_res2[res2])
                    elif res1 == res2:
                        intP = frequencies_res1[res1] * frequencies_res2[res2] * (1-CCplus_bias_I2_P_pcc) + frequencies_res1[res1] * (CCplus_bias_I2_P_pcc)
                        intP_bias = frequencies_res1[res1] * frequencies_res2[res2]
                    col_num = frequency_matrix_PCC[0].index(''.join(sorted(int)))
                    frequency_matrix_PCC[row_index][col_num] += intP
            elif int_code == "I3_P":
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = (frequencies_res1[res1]*frequencies_res2[res2]) * (1-CCplus_bias_I3_P_pcc)
                    elif res1 == res2:
                        #print(res1, frequencies[res1])
                        intP = (frequencies_res1[res1] * frequencies_res2[res2] * (1-CCplus_bias_I3_P_pcc)) + (frequencies_res1[res1] * (CCplus_bias_I3_P_pcc))
                        intP_bias = frequencies_res1[res1] * frequencies_res2[res2]
                        #print(intP, intP_bias)
                    col_num = frequency_matrix_PCC[0].index(''.join(sorted(int)))
                    frequency_matrix_PCC[row_index][col_num] += intP
            else:
                for int in interactions:
                    res1 = int[0]
                    res2 = int[1]
                    if res1 != res2:
                        intP = frequencies_res1[res1]*frequencies_res2[res2]
                    elif res1 == res2:
                        intP = frequencies_res1[res1] * frequencies_res2[res2]
                    if int_code == "I4" or int_code == "I1":
                        col_num = frequency_matrix_PCC[0].index(''.join(sorted(int)))
                    else:
                        col_num = frequency_matrix_PCC[0].index(int)
                    frequency_matrix_PCC[row_index][col_num] += intP

#saves expected probability matrix
frequency_df_PCC = pd.DataFrame(frequency_matrix_PCC[1:], columns = frequency_matrix_PCC[0])
frequency_df_APCC = pd.DataFrame(frequency_matrix_APCC[1:], columns = frequency_matrix_APCC[0])
frequency_df_PCC.to_csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/Expected_interaction_matrix_PCC_24_28.csv", index=False)
frequency_df_APCC.to_csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/Expected_interaction_matrix_APCC_24_28.csv", index=False)