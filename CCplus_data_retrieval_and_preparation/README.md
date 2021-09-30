Protocol CCplus
Step 1: Retrieve data from CCplus database (reduncancy : redundant, orientation: antiparallel/parallel, chains: any, a-helices: 2, Partnering: homodimer, repeats: canonical, Length > 11)
Then format coiled coils as <pdbID>,<coiledcoilID>,<start>,<finish>,<sequence>,<register> (1a59,0,95,106,DVARTAVSVLGA,defgabcdefga)
Step 2: put CCs in database and get information from the pdb database about the proteins (CCplus prepare.py) for CC_plus in SAMcc/CCplus *.txt file downloaded from the CC+ database as discribed above
Step 3: a. save unique bundles in fasta file by running the  the script (Unique_CC_fasta_for_cdhit.py)
	b. remove redundant AA-seqs using cd-hit at 50% threshold (http://weizhong-lab.ucsd.edu/cdhit_suite/cgi-bin/index.cgi?cmd=cd-hit)
		as bandwith of alignment the shortest length - 1 was taken
		as length of sequences to skip the shortest length - 1 was taken
		saved as text file under cdhit_output file under Input/Make_df_from_cdhit_fasta/cdhit_fastas using the same name as the output_file of Unique_CC_fasta_for_cdhit.py but replace last part (fasta_for_cdhit) with cdhit_out 
	c. execute the script (Make_df_from_cdhit_fasta.py) to create a cdhit_df safed in cdhit_dfs folder 
