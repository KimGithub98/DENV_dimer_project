# Script title: Interaction matrices GA output
# Version title: N/A
# Short application: creates interaction matrices of the GA output to identify which interactions were enriched by the GA
# Input files: csv file of the final GA population
# Output files: csv files with interaction matrices
# How to use: Run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 30/09/2021

import pandas as pd
import itertools
import math
import pathlib

hep_pos = "defgabcdefgabcdefgabcdefga"

def permutations_with_replacement(n, m):
    for i in itertools.product(m, repeat=n):
        yield i

def listToString(s):
    str1 = ""
    return (str1.join(s))

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

def empty_interaction_count_matrix():
    ## make empty interaction count matrix
    aalist = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    interactions = permutations_with_replacement(n = 2, m = aalist)
    interactions = [''.join(i) for i in interactions]
    interactions.insert(0, "pos")
    int_count_matrix = [interactions]
    int_list = ["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"]
    for int in int_list:
        listadd = [int, *itertools.repeat(0, 400)]
        int_count_matrix.append(listadd)
    return int_count_matrix

def get_int_count_matrix(table_down_up, int_count_matrix):
    for i in range(1, 4):
        for i_int in range(0, 4):
            interaction = listToString(sorted((table_down_up[0][i][i_int], table_down_up[1][i][i_int])))
            int_code = table_down_up[0][0][i_int]
            if interaction in int_count_matrix[0]:
                col = int_count_matrix[0].index(interaction)
                row = ["I2", "I3", "I1", "I4"].index(int_code) + 1
                int_count_matrix[row][col] = int_count_matrix[row][col] + 1
    return int_count_matrix

def get_up_down_helix(aa_code, heptad, int_count_matrix, orientation):
    pos_down = 0
    row = 1
    table_down = [["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"]]
    startpos = int("abcdefg".index(heptad[0]))
    rownumber = math.ceil((len(aa_code)+startpos)/7)
    for i in range(0, rownumber):
        listadd = [*itertools.repeat("0", 8)]
        table_down.append(listadd)

    table_up = [["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"]]
    startpos = int("gfedcab".index(heptad[::-1][0]))
    rownumber = math.ceil((len(aa_code)+startpos)/7)
    for i in range(0, rownumber):
        listadd = [*itertools.repeat("0", 8)]
        table_up.append(listadd)

    for aa in aa_code:
        hp = heptad[pos_down]
        down_dict = {"a": ['I2_P', 'I2_A'], "b": "I0", "c": "I0", "d": ["I3_P", "I3_A"], "e": ["I1", "I5"], "f": "I0", "g": ["I4", "I6"]}
        int_code = down_dict[hp]
        if int_code != "I0":
            if orientation == "PCC":
                up_dict = {"I2_P": "a", "I3_A": "a", "I6": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I5": "d", "I4": "e", "I0": "f", "I1": "g"} #parallel
                heptad_OR = heptad
                aa_code_OR = aa_code
            elif orientation == "APCC":
                up_dict = {"I2_P": "a", "I3_A": "a", "I5": "a", "I0": "b", "I0": "c", "I3_P": "d", "I2_A": "d", "I6": "d", "I1": "e", "I0": "f", "I4": "g"} #anti-parallel
                heptad_OR = heptad[::-1]
                aa_code_OR = aa_code[::-1]
            for int_c in int_code:
                up_hp = list(up_dict.values())[list(up_dict.keys()).index(int_c)]
                positions_up = [i for i, x in enumerate(heptad_OR) if x == up_hp]
                dis_small_4 = [abs(x - pos_down) < 4 for x in positions_up]
                if any(dis_small_4):
                    index = dis_small_4.index(True)
                    pos_up = positions_up[index]
                    column = table_down[0].index(int_c)
                    if int_c == "I2_P" or int_c == "I3_P":
                        interaction = ''.join(sorted(listToString((aa, aa_code_OR[pos_up]))))
                    elif int_c == "I1" or int_c == "I4":
                        interaction = ''.join(sorted(listToString((aa, aa_code_OR[pos_up]))))
                    else:
                        interaction = listToString((aa, aa_code_OR[pos_up]))
                    #if interaction == "EE":
                        #print(aa_code, pos_down)
                    if interaction in int_count_matrix[0]:
                        col = int_count_matrix[0].index(interaction)
                        row = ["I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6"].index(int_c) + 1
                        int_count_matrix[row][col] = int_count_matrix[row][col] + 1
        pos_down = pos_down + 1
    return int_count_matrix

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

CC_df = pd.read_csv('E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/run_final.csv', header=0)
int_count_matrix_APCC = empty_interaction_count_matrix()
int_count_matrix_PCC = empty_interaction_count_matrix()
i = 0
for CC in CC_df.aa_seq:
    i += 1
    if is_period(hep_pos):
        for orientation in ["PCC", "APCC"]:
            if orientation == "APCC":
                int_count_matrix_APCC = get_up_down_helix(aa_code=CC, heptad=hep_pos, int_count_matrix=int_count_matrix_APCC, orientation= "APCC")
            elif orientation == "PCC":
                int_count_matrix_PCC = get_up_down_helix(aa_code=CC, heptad=hep_pos, int_count_matrix=int_count_matrix_PCC, orientation="PCC")

for matrix in [int_count_matrix_PCC, int_count_matrix_APCC]:
    name_matrix = namestr(matrix, globals())[0]
    df = pd.DataFrame(matrix[1:], columns=matrix[0])
    del df['pos']
    df["Sum"] = df.sum(axis=1)
    print(df)
    df_P = df.div(df["Sum"], axis=0)
    print(df_P)
    pathlib.Path("E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/Interaction matrix/").mkdir(parents=True, exist_ok=True)
    df.to_csv("E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/Interaction matrix/_as" + name_matrix[16:] + ".csv", index=False)
    pathlib.Path("E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/Interaction matrix/P/" ).mkdir(parents=True, exist_ok=True)
    df_P.to_csv("E:/DENV_dimer_project/Output/Genetic Algorithm/Final_run/Interaction matrix/P/_as" + name_matrix[16:]+ "_P" + ".csv", index=False)
