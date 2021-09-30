import pandas as pd
import pypdb
import re


# Python program to convert a list
def listToString(s):
    str1 = " "
    return (str1.join(s))


# import files
pdb_scan = pd.read_csv("e:/SAMcc/pdb_scan_topology_segments_20200717.csv", header=0)
homodimers = pd.read_csv("e:/SAMcc/output/antiCC_homodimers_dictR.csv", header=None)
# print(homodimers)
print(pdb_scan)

# make dataframe by integrating all information of homodimers

homodimerdf = pd.DataFrame(columns=["pdbID", "bundleID", "Molecule", "gene", "organism", "CC_length", "from_p", "to_p", "chain1", "chain2", "pub_Title"])

for bundle in homodimers[0]:
    pdb = bundle[0:4]
    bundleID = bundle[:-4]
    print(bundleID)
    pdbrow = pdb_scan.loc[pdb_scan['bundleid'] == bundleID]
    length = pdbrow.length
    from_p = pdbrow.from_P
    to_p = pdbrow.to_P
    chain1 = homodimers[1][homodimers[0] == bundle]
    chain2 = homodimers[2][homodimers[0] == bundle]
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
    homodimerlist = [pdb, bundleID, molecule, gene, organism, length.values[0], from_p.values[0], to_p.values[0],
                     chain1.values[0], chain2.values[0], pub_title]
    print(homodimerlist)
    df_length = len(homodimerdf)
    homodimerdf.loc[df_length] = homodimerlist

print(homodimerdf)
homodimerdf.to_csv("e:/SAMcc/output/antiCC_homodimers_dfR.csv", index=False)
