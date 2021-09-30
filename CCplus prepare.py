# Script title: CCplus prepare
# Version title: N/A
# Short application: prepares retrieved data from the CCplus database and creates a workable .CSV file with coiled coil information and pdb information.
# Input files: txt file retrieved from CCplus database with format: <pdbID>,<coiledcoilID>,<start>,<finish>,<sequence>,<register> (example files: ~/DENV_dimer_project/Input/CCplus_database/Used_and_examples)
# Output files: csv file with format: <pdbID>,<bundleID>,<Molecule>,<gene>,<CC_length>,<from_p>,<to_p>,<chain1>,<chain2>,<poschain1>,<poschain2>.<pud_Title>
# How to use: Specify path from (E:/DENV_dimer_project/Input/CCplus_database) to the input file under user_input. Check in text file if the lines start with Tabs, if this is the case: uncommend pdb = pdb[6:10] line 27

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 16/07/21

user_input = "Used_and_examples/CCplus_PCC" #insert string of name and location of the txt input file in ~/DENV_dimer_project/Input/CCplus_database folder. Example: "Used_and_examples/CCplus_APCC"

import pandas as pd
import pypdb
import re

CCplus_APCC = pd.read_csv("E:/DENV_dimer_project/Input/CCplus_database/" + user_input + ".txt", header=None)
def listToString(s):
    str1 = " "
    return (str1.join(s))

CCplus_APCCdf = pd.DataFrame(columns=["pdbID", "bundleID", "Molecule", "gene", "organism", "CC_length", "from_p", "to_p", "chain1", "chain2", "poschain1", "poschain2", "pub_Title"])

for i in range(0, len(CCplus_APCC)):
    pdb = CCplus_APCC[0][i]
    pdb = pdb[6:10] #uncommend if tabs in file
    bundle = CCplus_APCC[1][i]
    bundleID = str(pdb) + "_" + str(bundle)
    end = CCplus_APCC[3][i]
    start = CCplus_APCC[2][i]
    length = int(end)-int(start)
    chain1 = CCplus_APCC[4][i]
    pos_chain1 = CCplus_APCC[5][i]
    pdb_file = pypdb.get_pdb_file(pdb, filetype='pdb', compression=False)
    pd_split = pdb_file.splitlines()
    terms = ["TITLE", "COMPND   2", "SOURCE   2 ORGANISM_SCIENTIFIC: ", "SOURCE   5 GENE: "]
    title = [s for s in pd_split if terms[0] in s]
    pub_title = re.sub(r"TITLE\s{5}|TITLE\s{4}\d\s|\s{2,}", "", listToString(title))
    molecule = [s for s in pd_split if terms[1] in s]
    molecule = re.sub(r"COMPND   2 MOLECULE: |;\s*", "", listToString(molecule))
    organism = [s for s in pd_split if terms[2] in s]
    organism = re.sub(r"SOURCE   2 ORGANISM_SCIENTIFIC: |;\s*", "", listToString(organism))
    gene = [s for s in pd_split if terms[3] in s]
    gene = re.sub(r"SOURCE   5 GENE: |;\s*", "", listToString(gene))
    list = [pdb, bundleID, molecule, gene, organism, length, "NA", "NA", chain1, "NA", pos_chain1, "NA", pub_title]
    df_length = len(CCplus_APCCdf)
    CCplus_APCCdf.loc[df_length] = list

print(CCplus_APCCdf)
CCplus_APCCdf.to_csv("E:/DENV_dimer_project/Output/CCplus_database/" + user_input + ".csv", index=False)
