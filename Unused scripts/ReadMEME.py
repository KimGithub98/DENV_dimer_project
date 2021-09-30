import re
import pandas as pd
import os

#os.chdir("/home/nmrbox/0013/kwijnant/Downloads/meme-5.3.3")

blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}
CC_df_in = pd.read_csv("e:/SAMcc/output/cdhit_dfs/cdhit_length14_18_APCC_HOMO_df.csv", header=0)
CC_df_in2 = pd.read_csv("e:/SAMcc/output/cdhit_dfs/cdhit_CCplus_length14_18_APCC_HOMO_df.csv", header=0)
CC_df_in = pd.concat([CC_df_in, CC_df_in2])
#CC_df_in = pd.read_csv("cdhit_length14_18_APCC_HOMO_df.csv", header=0)

def get_motif_scores(sites, CC_df, motif_inf):
    sites = sites.split("\n")
    hep_starts = {"a": 0, "b": 0, "c": 0, "d": 0, "e": 0, "f": 0, "g": 0}
    print(sites[1])
    match = re.search("sites =[\s|\d]{4}", motif_inf)
    total = int(re.search("\d+", match[0])[0])
    for line in sites[5:(total + 5)]:
        line = re.split("(\d)\s+", line)
        line[0:2] = [''.join(line[0:2])]
        line[3:5] = [''.join(line[3:5])]
        line.pop(1)
        start = int(line[1])
        bundle_id = line[0]
        hep_pos = CC_df.poschain1[CC_df.bundleID == bundle_id]
        hep_pos = hep_pos.values[0]
        hep_start = hep_pos[(start-1)]
        hep_starts[hep_start] += 1
    print(sorted(hep_starts.values(),reverse=True)[1])
    sites_score = (max(hep_starts.values())-sorted(hep_starts.values(),reverse=True)[1])/sum(hep_starts.values())
    print(hep_starts)
    print("score is: " + str(sites_score))
    pred_hep_pos = list(hep_starts.keys())[list(hep_starts.values()).index(max(hep_starts.values()))] #return
    return (hep_starts, pred_hep_pos, sites_score)

def get_aa_scores(description, start_hep_pos, motif_inf):
    description = description.split("\n")
    i = 0
    for line in description:
        if "Multilevel" in line:
            line = re.split("\s{11}", line)
            Multilevel = line[1]
        elif 3 <= i <= 22:
            line = line[18::]
            line = re.split("(\D)\s{2}", line)
            aa = line[1]
            count = line[2]
            count = re.sub(":", "0", count)
            PSPM_add = pd.DataFrame({aa: list(count)})
            if i == 3:
                PSPM = pd.DataFrame.from_dict(PSPM_add)
            else:
                PSPM = PSPM.join(PSPM_add)
        i += 1
    pos = 0
    score_list = []
    dict_hep_pos = {"a": 0, "b": 1, "c": 2, "d": 3, "e": 4, "f": 5, "g": 6}
    hep_pos = dict_hep_pos[start_hep_pos]
    for aa in Multilevel:
        if hep_pos in [0,3,5,6]:
            pos_PSPM = PSPM.iloc[[pos],]
            match = re.search("sites =[\s|\d]{4}", motif_inf)
            total = int(re.search("\d+", match[0])[0])
            aa_score = 0
            for substitution in pos_PSPM:
                count = pos_PSPM[substitution].item()
                if count == "a": count = total
                if int(count) > 0:
                    if (aa, substitution) in blosum62.keys():
                        score = blosum62[(aa, substitution)]
                    else:
                        score = blosum62[(substitution, aa)]
                    score_new = (int(count)/total)*score
                    #print("count: " + str(count) + " total: " + str(total) + " new score: "  + str(score_new))
                    aa_score += score_new
            score_list.append(aa_score)
            #print("score: " + str(aa_score))
        pos += 1
        hep_pos += 1
        if hep_pos == 7:
            hep_pos = 0
    print(score_list)
    final_score = sum(score_list)/len(score_list)
    print(final_score)
    return final_score

#"meme_out/meme.txt"
#"e:/SAMcc/meme.txt"
with open("e:/SAMcc/meme4.txt") as fo:
    print(fo)
    op = ""
    start = 0
    i = 1
    ignore_next_line = False
    cntr1 = 0
    cntr2 = 0
    Description = ""
    sites = ""
    final_Description = ""
    final_sites = ""
    sites_information = False
    Des_information = False
    for x in fo.read().split("\n"):
        if ignore_next_line:
            x = ""
            ignore_next_line = False
        if x == "********************************************************************************":
            motif = "new"
            if start == 1:
                # print(op)
                # print(final_sites)
                # print(repr(final_Description))
                motif_score = get_motif_scores(sites=final_sites, CC_df=CC_df_in, motif_inf= motif_header)
                if motif_score[2] >= 0.5:
                    aa_score = get_aa_scores(description=final_Description, start_hep_pos= motif_score[1], motif_inf= motif_header)
                    print("the final motif score = " + str(aa_score*2/(1/motif_score[2])**2))
                op = ""
                i = i + 1
            else:
                start = 1
                op = ""
                i = 1
        else:
            if motif == "new" and "MOTIF " in x:
                motif_header = x
                op = op + "\n" + x
                motif = "started"
                ignore_next_line = True
            elif motif == "started":
                op = op + "\n" + x
                if "Description" in x or Des_information is True:
                    Description = Description + "\n" + x
                    Des_information = True
                    if x == "--------------------------------------------------------------------------------":
                        if cntr1 == 0:
                            # ignores first dashed line
                            cntr1 += 1
                        else:
                            # end special information
                            # print(Description)
                            final_Description = Description
                            cntr1 = 0
                            Description = ""
                            Des_information = False
                elif "sites sorted by position p-value" in x or sites_information is True:
                    sites = sites + "\n" + x
                    sites_information = True
                    if x == "--------------------------------------------------------------------------------":
                        if cntr2 == 0:
                            # ignores first dashed line
                            cntr2 += 1
                        else:
                            # end special information
                            # print(sites)
                            final_sites = sites
                            cntr2 = 0
                            sites = ""
                            sites_information = False
                else:
                    pass
            else:
                start = 0
fo.close()

import sys
print(sys.version)