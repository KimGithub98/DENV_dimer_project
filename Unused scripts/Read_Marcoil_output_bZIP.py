## creates CSV in >1yke_1_2,abcdefgabcdefga,IMLLEEQLEYKRGEI, lowest Marcoil score, Max Marcoil score, Threshold 10 CC length, Threshold 50 CC length format
import re
import pandas as pd
import numpy as np
import os

run_name = "bZIP"
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

def true_in_a_row(bool_list):
    i = 0
    i_max = 0
    for bool in bool_list:
        if bool == True:
            i +=1
            if i > i_max:
                i_max = i
        elif bool == False:
            i = 0
    return i_max

for run in [run_name]:
    for file in os.listdir('e:/SAMcc/output/genetic_algorithm_output/Marcoil_out'):
        if run in file:
            mar_out = open('e:/SAMcc/output/genetic_algorithm_output/Marcoil_out/' + file)
    cntr = 0
    new_CC = False
    results = False
    register = ""
    aa_seq = ""
    prob_list = []
    Matrix = [["name", "register", "AA sequence", "min score (MARCOIL)", "max score (MARCOIL)", "average score (MARCOIL)", "Length CC Threshold 0", "Length CC Threshold 10", "Length CC Threshold 50"]]
    for x in mar_out.read().split("\n"):
        if "****" in x:
            # save CC
            prob_array = np.array(prob_list)
            Lenght_T_0 = true_in_a_row(prob_array > -1)
            Lenght_T_10 = true_in_a_row(prob_array > 10)
            Lenght_T_50 = true_in_a_row(prob_array > 50)
            list = [CC_name, register, aa_seq, min(prob_list), max(prob_list), sum(prob_list)/len(prob_list), Lenght_T_0, Lenght_T_10, Lenght_T_50]
            Matrix.append(list)
            # empty all parameters
            aa_seq = ""
            register = ""
            prob_list = []
            new_CC = False
            results = False
            cntr = 0
        elif ">" in x:
            new_CC = True
            cntr += 1
            CC_name = re.search(r'>\S*\s', x)[0]
        elif new_CC is True and cntr == 1:
            x = x.replace(" ", "")
            #aa_seq = x.replace("*", "")
            cntr += 1
        elif results is True and cntr > 1:
            x_search = re.search(r'\s*(\d*)\s(\S)\s*(.*)\s(\S)', x)
            if x_search and float(x_search.group(3)) > 1:
                aa_num = x_search.group(1)
                aa_seq = aa_seq + x_search.group(2)
                prob_list.append(float(x_search.group(3)))
                register = register + x_search.group(4)
        elif "cc-probability in percent and best heptad phase" in x:
            results = True
            # ignore this line but save next lines
    df = pd.DataFrame(Matrix[1:], columns= Matrix[0])
    pd.set_option('display.max_columns', None)  # or 1000
    pd.set_option('display.max_rows', 10)  # or 1000
    pd.set_option('display.max_colwidth', None)  # or 199
    df.to_csv("e:/SAMcc/bZIP affinity test/MARCOIL.csv", index=False)
    #print(df.head)

minor_helical_phase = {"a": "197", "b": "300", "c": "41", "d": "146", "e": "249", "f": "351", "g": "95"}
for name in df["name"]:
    register = df["register"][df["name"] == name]
    sequence = df["AA sequence"][df["name"] == name]
    print(len(register.values[0]), len(sequence.values[0]))
    phase = minor_helical_phase[register.values[0][0]]
    #print(phase)
    if is_period(register.values[0]):
        parameters_file = open("e:/SAMcc/bZIP affinity test/parameters_structure/" + name[1:len(name) - 1] + ".txt",
                               "w")
        sequence_file = open("e:/SAMcc/bZIP affinity test/sequence_files_half/" + name[1:len(name) - 1] + ".txt",
                               "w")
        register_file = open("e:/SAMcc/bZIP affinity test/register_files_half/" + name[1:len(name) - 1] + ".txt",
                             "w")
        parameters_text = """4helix # base name of output structures\n""" +\
                          str(len(register.values[0])) + """ # number of residues
    1.82 1.82 1 # twist range (low, high, number of samples)
    4.9 4.9 1 # supercoil radius range (low,high, number of samples)
    2 # number of helices
    2 # number of helices in asymetric unit (an entry for each helix is required below)
    1 1 # orientation of each helix in asymetric unit\n""" + \
                          phase + " " + phase + " 1 # phase of first helix (low, high, number of samples)\n" + phase + " " + phase + """ 1 # phase of second helix (low, high, number of samples)
    0 0 1 # z-offset of first helix (low, high, number of sampels)
    0 0 1 # z-offset of second helix (low, high, number of sampels)
    A B # chain names
    0 1 """
        #print(parameters_text)
        parameters_file.write(parameters_text)
        sequence_file.write(sequence.values[0])
        #sequence_file.write(sequence.values[0])
        register_file.write(register.values[0])
        parameters_file.close()
        sequence_file.close()
        register_file.close()

