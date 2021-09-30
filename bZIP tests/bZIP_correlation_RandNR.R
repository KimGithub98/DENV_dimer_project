# Script title: correlation bZIPs and (non) relative scores
# Version title: N/A
# Short application: calculates correlation between affinity and scores
# Input files: .csv of the (non)relative scores of the bZIPs and .csv with kd values of the bZIPs
# Output files: csv file containing the correlation of the scores and bZIPs
# How to use: specify if you want to calculate correlations of relative or non-relative scores at score = under user input and run script.

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 29/09/2021

library("ggpubr")
library("ggExtra")

#USER_INPUT
score_name = "Relative" #choose between Relative/Non_Relative

#fetch input files
if (score_name == "Relative"){
  bZIP_relative_scores <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/bZIP_Relative_Scores.csv")
} else if (score_name == "Non_Relative"){
  bZIP_relative_scores <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/bZIP_Non_Relative_Scores.csv")
}

bZIP_affinity <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/KD_and_sequence_files/kd.csv", header=TRUE)
colnames(bZIP_affinity) <- c("pdb", "pdb2", "G", "Length")

#creates a dataframe of the scores and how they are correlated with affinity
cor_df <- data.frame(matrix(ncol = 3, nrow = 0))
i = 0
mean_length <- mean(bZIP_affinity$Length)
for (plot_score in colnames(bZIP_relative_scores)[-c(1:6)]){
  i = i + 1
  plot_df = data.frame(0, 0)
  colnames(plot_df) = c("G", "score")
  # looks at scores and Kd values of per pdb 
  for (pdb in bZIP_relative_scores[,1]){
    G_free = bZIP_affinity[bZIP_affinity$pdb ==pdb, "G"]
    Score = bZIP_relative_scores[bZIP_relative_scores$pdb==pdb, plot_score]
    length_comp = mean(bZIP_relative_scores$num_int_APCC)/bZIP_relative_scores[bZIP_relative_scores$pdb==pdb, "num_int_APCC"] #used to compensate for the increase in Kd because of more interactions/length
    if (!is.null(G_free)){
      plot_df = rbind(plot_df, data.frame("G" = (((293.15*8.31446261815324)*log(as.numeric(G_free)/1e9)/4186.8))/-1, "score"= abs(as.numeric(Score))*length_comp)) #calculates G from Kd
    }
  }
  colnames(plot_df) = c("G", "score")
  plot_df = plot_df[-1,]
  test_prod <- cor.test(plot_df$G, plot_df$score) #calculates correlation
  cor_df = rbind(cor_df, c(plot_score, test_prod$estimate, test_prod$p.value))
}

#saves dataframe
colnames(cor_df) = c("score_component", "R-value", "p-value")
write.csv(cor_df, paste("E:/DENV_dimer_project/Output/bZIP tests/", score_name, "_score_bZIP_correlation.csv", sep = ""), row.names = FALSE)

