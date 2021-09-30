# Script title: Expected_interaction_matrix
# Version title: interaction_bias
# Short application: creates a dataframe that describes per length (24-28) and begin position (a-g) CC what the symmetry bias is for I2_P, I3_P for PCC and I1 and I4 for APCC. (1 is 100% symmetry and 0 is 0% symmetry)
# Input files: NA
# Output files: a csv file with the length and startposition (24a, 24b etc.) as rows and interaction types as columns
# How to use: Change length_pos to the length range you want and run

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 23/09/2021

import itertools
import pandas as pd

length_pos = [24, 25, 26, 27, 28] #change to different length range if wanted

## calculate interaction bias for every CC length and start heptad position
start_pos_pos = {"a":0, "b":1, "c":2, "d":3, "e":4, "f":5, "g":6}
heplist = ["a1", "b1", "c1", "d1", "e1", "f1", "g1", "a2", "b2", "c2", "d2", "e2", "f2", "g2", "a3", "b3", "c3", "d3", "e3", "f3", "g3", "a4", "b4", "c4", "d4", "e4", "f4", "g4", "a5", "b5", "c5", "d5", "e5", "f5", "g5"]
table_out = [["rownames", "I1", "I2_P_pcc", "I3_P_pcc", "I2_P_apcc", "I3_P_apcc", "I4"]]
table_out2 = [["rownames", "I1", "I2_P_pcc", "I3_P_pcc", "I2_P_apcc", "I3_P_apcc", "I4"]]
i2 = 1

for length in length_pos: #for every length
    for start_pos_key in start_pos_pos.keys(): #for every starting position of the heptad register
        start_pos = start_pos_pos[start_pos_key]
        heptad_sequence = []
        #creates heptad register of the CC of a length with a certain starting heptad position
        for i in range(0, length):
            pos = i + start_pos
            heptad_sequence.append(heplist[pos])

        #creates tables that describe the bias
        pos_down = 0
        itself = 0
        not_itself = 0
        table_out.append([str(length) + start_pos_key, 0, 0, 0, 0, 0, 0]) #empty count matrix for the  times a position binds to itself (symmetrically)
        table_out2.append([str(length) + start_pos_key, 0, 0, 0, 0, 0, 0]) #empty count matrix for the times a position binds to a other residue (non-symmetrically)
        for aa in heptad_sequence:
            hp = aa[0]
            down_dict = {"a": ['I2_P', 'I2_A'], "b": "I0", "c": "I0", "d": ["I3_P", "I3_A"], "e": ["I1", "I5"],
                         "f": "I0", "g": ["I4", "I6"]}
            int_code = down_dict[hp]
            #only takes into accout I2_P I3_P for PCC and I1 and I4 for APCC as these can have bias, the others can't as they bind to another heptad position
            if int_code != "I0":
                for int_c in int_code:
                    if int_c != "I5" or int_c != "I6" or int_c != "I2_A" or int_c != "I3_A":
                        for orientation in ["APCC", "PCC"]: #for each orientation
                            heptad_OR = "h" #non-existing position to avoid double counting
                            # the if statements below create a dictionary and heptad register based on whether its a APCC or PCC
                            if orientation == "APCC" and (int_c == "I2_P" or int_c == "I3_P"):
                                ## APCC
                                up_dict = {"I2_P": "a", "I3_A": "a", "I5": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I6": "d", "I1": "e", "I0": "f", "I4": "g"} #anti-parallel
                                heptad_OR = heptad_sequence[::-1] #referces heptad register as this is a APCC
                                int_c_o = int_c + "_apcc"
                            if orientation == "PCC" and (int_c == "I2_P" or int_c == "I3_P"):
                                ## PCC
                                up_dict = {"I2_P": "a", "I3_A": "a", "I6": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I5": "d", "I4": "e", "I0": "f", "I1": "g"} #parallel
                                heptad_OR = heptad_sequence
                                int_c_o = int_c + "_pcc"
                            if orientation == "APCC" and (int_c == "I1" or int_c == "I4"):
                                ## APCC
                                up_dict = {"I2_P": "a", "I3_A": "a", "I5": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I6": "d", "I1": "e", "I0": "f", "I4": "g"} #anti-parallel
                                heptad_OR = heptad_sequence[::-1] #referces heptad register as this is a APCC
                                int_c_o = int_c
                            #checks if position binds symmetricly to itself in a homodimer
                            up_hp = list(up_dict.values())[list(up_dict.keys()).index(int_c)] #gets heptad position of binding partner
                            #looks at which heptad position is closest to find the binding partner
                            positions_up = [i for i, x in enumerate(heptad_OR) if x[0] == up_hp]
                            dis_small_4 = [abs(x - pos_down) < 4 for x in positions_up]
                            #assesses if its binds to itself and counts the times this happens for a interaction
                            column_num = {"I1":1, "I2_P_pcc":2, "I3_P_pcc":3,"I2_P_apcc":4, "I3_P_apcc":5, "I4":6}
                            if any(dis_small_4):
                                up_hp_code = list(itertools.compress(positions_up, dis_small_4))
                                up_hp_code = up_hp_code[0]
                                if aa == heptad_OR[up_hp_code]:
                                    #they bind to themselfs
                                    table_out[i2][column_num[int_c_o]] += 1
                                else:
                                    #they do not bind to itself
                                    table_out2[i2][column_num[int_c_o]] += 1
            pos_down = pos_down + 1
        i2 += 1


df = pd.DataFrame(table_out)
df.columns = df.iloc[0]
df = df[1:]
df2 = pd.DataFrame(table_out2)
df2.columns = df2.iloc[0]
df2 = df2[1:]

#calculates the ratio between bind with itself and bind with not-itself
df3 = df[df.columns[1:]]/(df[df.columns[1:]]+df2[df2.columns[1:]])
df3.insert(0, "Coiled_coil", df[df.columns[0]])

##save interaction bias file
df3.to_csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/interaction_bias_24_28.csv", index=False)




