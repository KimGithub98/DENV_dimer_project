# import files
pdb_scan = open("e:/SAMcc/pdb_scan_nr_topology_segments_20200717.csv", "r") ##pdb_scan_nr_topology_segments_20200717.csv

#get AntiCCs
all_antiCC = []
for bundle in pdb_scan:
    bundle_split = bundle.split(',')
    top = bundle_split[6] #number of chains and bundle type 2_0 is antiparallel CC dimer
    from_p = bundle_split[12]
    to_p = bundle_split[13]
    is_antiCC = top == "2_0"
    if is_antiCC:
        all_antiCC.append(bundle_split[1])

#print(len(all_antiCC)) # give 5572 antiCCs
#print(all_antiCC)
print('2lw9_0' in all_antiCC)

#get homodimer/heterodimer
import os
import pandas as pd
from itertools import compress

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'MSE': 'M', 'UNK': 'X'}

def convert(s):
    str1 = ""
    return (str1.join(s))

homodimers = {}
csv = ".csv"
for bundles in all_antiCC:
    folder = bundles[1:3]
    filename = bundles + csv
    try:
        with open(os.path.join("e:/SAMcc/pdb_scan_csv_nonredundant", folder, filename)) as f:
            csv_file = pd.read_csv(f)
            #print(csv_file)
            chain1 = csv_file["chain"] == 0
            chain2 = csv_file["chain"] == 1
            resname = csv_file["res_name"]
            name_chain1 = list(compress(resname, chain1))
            name_chain2 = list(compress(resname, chain2))
            name_chain2.reverse()
            str1 = convert(name_chain1)
            str2 = convert(name_chain2[1:-1])
            if str2 in str1:
                #print("MATCH", filename, name_chain1, name_chain2)
                seq_chain1 = ""
                seq_chain2 = ""
                for item in name_chain1:
                    aa = d[item]
                    seq_chain1 = seq_chain1 + aa
                for item in name_chain2:
                    aa = d[item]
                    seq_chain2 = seq_chain2 + aa
                print(seq_chain2)
                hep_pos = csv_file["positions_samcc"]
                hep_pos_chain1 = convert(list(compress(hep_pos, chain1)))
                hep_pos_chain2 = convert(list(compress(hep_pos, chain2)))
                homodimers[filename] = [seq_chain1, seq_chain2, hep_pos_chain1, hep_pos_chain2]



    except IOError:
        print(filename, " not found")

#print(homodimers)
#print(homodimers.keys())

import pandas as pd

(pd.DataFrame.from_dict(data=homodimers, orient='index')
   .to_csv("e:/SAMcc/output/antiCC_homodimers_dictNR2 .csv", header=False))




