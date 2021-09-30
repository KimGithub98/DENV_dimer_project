## creates CSV in >1yke_1_2,abcdefgabcdefga,IMLLEEQLEYKRGEI, lowest Marcoil score, Max Marcoil score, Threshold 10 CC length, Threshold 50 CC length format
import re
import pandas as pd
import numpy as np
import os

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

for run in ["run0"]:
    gen_al_out = run + ".csv"
    for file in os.listdir('e:/SAMcc/output/genetic_algorithm_output'):
        if gen_al_out in file:
            print(file)
            gen_al_out_df = pd.read_csv('e:/SAMcc/output/genetic_algorithm_output/' + file, header = 0)
    for file in os.listdir('e:/SAMcc/output/genetic_algorithm_output/Marcoil_out'):
        if run in file:
            mar_out = open('e:/SAMcc/output/genetic_algorithm_output/Marcoil_out/' + file)
    cntr = 0
    new_CC = False
    results = False
    register = ""
    prob_list = []
    Matrix = [["name", "register", "AA sequence", "APCC_score", "CC_score", "Total_score", "min score (MARCOIL)", "max score (MARCOIL)", "average score (MARCOIL)", "Length CC Threshold 10", "Length CC Threshold 50"]]
    for x in mar_out.read().split("\n"):
        if "****" in x:
            # save CC
            prob_array = np.array(prob_list)
            Lenght_T_10 = true_in_a_row(prob_array > 10)
            Lenght_T_50 = true_in_a_row(prob_array > 50)
            list = [CC_name, register, aa_seq, "NA", "NA", "NA", min(prob_list), max(prob_list), sum(prob_list)/len(prob_list), Lenght_T_10, Lenght_T_50]
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
            aa_seq = x.replace("*", "")

            cntr += 1
        elif results is True and cntr > 1:
            x_search = re.search(r'\s*(\d*)\s(\S)\s*(.*)\s(\S)', x)
            if x_search:
                aa_num = x_search.group(1)
                prob_list.append(float(x_search.group(3)))
                register = register + x_search.group(4)
        elif "cc-probability in percent and best heptad phase" in x:
            results = True
            # ignore this line but save next lines
    df = pd.DataFrame(Matrix[1:], columns= Matrix[0])
    pd.set_option('display.max_columns', None)  # or 1000
    pd.set_option('display.max_rows', 10)  # or 1000
    pd.set_option('display.max_colwidth', None)  # or 199
    print(df)
    print(gen_al_out_df)

    for CC_name in df['name']:
        print(CC_name[1:5])
        print(gen_al_out_df['aa_seq'].values)
        if CC_name[1:5] in gen_al_out_df['aa_seq'].values:
            row1 = gen_al_out_df.loc[gen_al_out_df['aa_seq'] == CC_name[1:5]]
            df.loc[df['name'] == CC_name, "APCC_score"] = row1["APCC_score"].values[0]
            df.loc[df['name'] == CC_name, "CC_score"] = row1["CC_score"].values[0]
            df.loc[df['name'] == CC_name, "Total_score"] = row1["Total_score"].values[0]
            df.loc[df['name'] == CC_name, "register"] = row1["register"].values[0]
            row2 = df.loc[df['name'] == CC_name]
    print(df)
    df.to_csv("e:/SAMcc/output/genetic_algorithm_output/" + run + "_Marout_genout.csv", index=False)