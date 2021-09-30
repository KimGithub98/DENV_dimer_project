Protocol for general (CCplus) sequence motif generation and interaction matrices/plots
Sequence motifs
Step 1: Run sequence_logos_heptad.py and sequence_logos_full.py, this saves different sequence motifs and weight matrices in several layout folders in the folder sequence logos. This is done for all files in the cdhit_dfs folder.

Int_matrices
(choose all_pairs scripts for additional interaction pairs such as e-a' g-d')
Step 1: calculate interaction bias (interect with themselves at certain positions) to calculate expected interaction frequencies (Expected_interaction_matrix_interaction_bias.py). Edit length list for different length
Step 2: Make an interaction count matrix based on calculated interaction biases (Expected_interaction_matrix_make_probability_matrix.py)
Step 3: To make interaction profiles/matrices of a set of coiled coils in the cdhit_df folder run (Make_interaction_matrices.py)

Int_plots/csvs make (R scripts)
(choose all_pairs scripts for additional interaction pairs such as e-a' g-d')
Step 1: Run Interactions_volcanoplots(_all_pairs).R to save in E:samcc/output/Interaction volcano plots (this is done for one set of CCs)
Step 2: For rel probability matrices run Interactions_csvfile(_all_pairs).R, for non relative probability matrix run Non-relative_interaction_csvfile.R 
