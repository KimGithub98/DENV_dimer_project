# Script title: plot correlation bZIPs and final scores
# Version title: N/A
# Short application: calculate correlation between final scores and bZIPs and plot on a correlation plot
# Input files: .csv of the final scores of the bZIPs and .csv with kd values of the bZIPs
# Output files: .svg of the correlation between affinity, Marcoil and logicoil and the 3 elements of the total score and total score
# How to use: run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 29/09/2021

library("ggpubr")
library("ggExtra")

bZIP_final_scores <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/bZIP_final_scoring.csv")
bZIP_affinity <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/KD_and_sequence_files/kd.csv", header=TRUE)
colnames(bZIP_affinity) <- c("pdb", "pdb2", "G", "Length")

plot_list <- vector('list', 9)
cor_df <- data.frame(matrix(ncol = 3, nrow = 0))
cor_df_marcoil <- data.frame(matrix(ncol = 3, nrow = 0))
cor_df_logi <- data.frame(matrix(ncol = 6, nrow = 0))
i = 0
mean_length <- mean(bZIP_affinity$Length)

#colnames(bZIP_final_scores)[c(5,54,55,56)]
for (plot in c("G", "marcoil", "specificity")){
  for (plot_score in colnames(bZIP_final_scores)[c(4,5,6,7)]){
    ## for plot to compare LOGICOIL scores with APCC score
    plot_df = data.frame(0, 0, 0, 0)
    colnames(plot_df) = c("G", "score", "marcoil", "specificity")
    for (pdb in bZIP_final_scores[,1]){
      G_free = bZIP_affinity[bZIP_affinity$pdb ==pdb, "G"]
      Marcoil = MARCOIL_out[MARCOIL_out$name == paste(">",pdb, " ", sep = ""),"average.score..MARCOIL."]
      Logi_APCC = LOGI_MAR_bZIP[LOGI_MAR_bZIP$id == paste(">",pdb, " ", sep = ""), "Antiparallel.dimer"]
      Logi_PCC = LOGI_MAR_bZIP[LOGI_MAR_bZIP$id == paste(">",pdb, " ", sep = ""), "Parallel.Dimer"]
      Score = bZIP_final_scores[bZIP_final_scores$pdb==pdb, plot_score]
      length_comp = mean(bZIP_final_scores$num_int_APCC)/bZIP_final_scores[bZIP_final_scores$pdb==pdb, "num_int_APCC"]
      if (!is.null(G_free)){
        plot_df = rbind(plot_df, data.frame("G" = (((293.15*8.31446261815324)*log(as.numeric(G_free)/1e9)/4186.8))/-1, "score"= abs(as.numeric(Score))*length_comp, "marcoil" = Marcoil, "specificity" = as.numeric(Logi_PCC)-as.numeric(Logi_APCC)))
      }
    }
    colnames(plot_df) = c("G", "score", "marcoil", "specificity")
    plot_df = plot_df[-1,]
    
    test_prod <- cor.test(plot_df[,plot], plot_df$score)
    cor_num_prod <- paste("R = ", as.character(round(test_prod$estimate, digits = 4)))
    cor_df = rbind(cor_df, c(plot_score, test_prod$estimate, test_prod$p.value))
    cor_df_marcoil = rbind(cor_df_marcoil, c(plot_score, test_prod2$estimate, test_prod2$p.value))
    cor_df_logi = rbind(cor_df_logi, c(plot_score, test_prod3$estimate, test_prod3$p.value))
    
    i = i + 1
    p <- ggscatter(plot_df, x = plot, y = "score", 
                   size = 0.3, alpha = 0.5, 
                   add = "reg.line", conf.int = TRUE, cor.method = "pearson",
                   xlab = as.character(plot), ylab = as.character(plot_score)) + annotate(geom = "text", x = max(plot_df[,plot]), y = max(plot_df[,"score"]), label = cor_num_prod, hjust = 1, vjust = 1, size = 2)
    p <- p + font("xlab", size = 5) + font("ylab", size = 5) + font("xy.text", size = 7)
    p <- p + theme(aspect.ratio=1) 
    p <- p + theme_bw()
    p <- p + theme(panel.background = element_rect(fill = "grey95"),
                   text = element_text(family = "sans", colour = "#3E4648", size = 10), 
                   plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
                   aspect.ratio = 1/1.4,
                   plot.margin = unit(c(0,0,0,0), "cm"))
    plot_list[[i]] <- p
  }
}

colnames(cor_df) = c("score_component", "R-value", "p-value")
colnames(cor_df_marcoil) = c("score_component", "R-value", "p-value")
colnames(cor_df_logi) = c("score_component", "R-value (APCC)", "p-value (APCC)")
print(ggarrange(plotlist = plot_list, ncol = 4, nrow = 3))
ggsave(filename = "bZIP_final_scores_FINAL(total).svg", path = "E:/SAMcc/output/HADDOCK_CC_tests")
