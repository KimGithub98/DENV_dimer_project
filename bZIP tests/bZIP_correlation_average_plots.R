# Script title: scatter plot affinity correlation and average score
# Version title: N/A
# Short application: plots the Rvalue of affinity of a score against the average of that score in bZIPs
# Input files: csv of correlation and Relative scores
# Output files: .svg with the plot
# How to use: Change score_name in USER_INPUT to Relative or Non_Relative and run script

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 29/09/2021

library(ggplot2)

#USER_INPUT
score_name = "Non_Relative" #choose between Relative/Non_Relative

#fetch input files
if (score_name == "Relative"){
  bZIP_relative_scores <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/bZIP_Relative_Scores.csv")
  R_values <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/Relative_score_bZIP_correlation.csv")
} else if (score_name == "Non_Relative"){
  bZIP_relative_scores <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/bZIP_Non_Relative_Scores.csv")
  R_values <- read.csv("E:/DENV_dimer_project/Output/bZIP tests/Non_Relative_score_bZIP_correlation.csv")
}

#create a df for the plot
mean_scores <- data.frame(matrix(ncol = 2, nrow = 0))
for (score in colnames(bZIP_relative_scores)[7:38]){
  scores = bZIP_relative_scores[,score]
  average = mean(scores)
  range = c(min(scores), max(scores))
  mean_scores <- rbind(mean_scores, c(average, range))
}
mean_scores <- cbind(mean_scores, c("I1", "I2_P", "I2_A", "I3_P", "I3_A", "I4", "I5", "I6", "I1", "I2_P", "I2_A", "I3_P", "I3_A", "I4", "I5", "I6", "I1", "I2_P", "I2_A", "I3_P", "I3_A", "I4", "I5", "I6", "I1", "I2_P", "I2_A", "I3_P", "I3_A", "I4", "I5", "I6"))
mean_scores <- cbind(mean_scores, c("red","red","red","red","red","red","red","red", "green","green","green","green","green","green","green","green", "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow", "blue","blue","blue","blue","blue","blue","blue","blue"))
mean_scores <- cbind(mean_scores, R_values[1:32,2])
rownames(mean_scores) <- colnames(bZIP_relative_scores)[7:38]
colnames(mean_scores) <- c("average", "min", "max", "int", "outline", "R_value")
mean_groups <- data.frame(matrix(ncol = 2, nrow = 0))
for (group in unique(mean_scores$outline)){
  mean_average = mean(mean_scores$average[mean_scores$outline == group])
  mean_groups = rbind(mean_groups, c(mean_average,0))
}
mean_groups = cbind(mean_groups, c("I1", "I1", "I1", "I1"), c("red", "green", "yellow", "blue"))
colnames(mean_groups) = c("average", "R_value", "int", "outline")
rownames(mean_groups) = c("APCC1", "PCC1", "APCC2", "PCC2")

#create a plot
pointcols<- c("chartreuse3", "red4", "brown2", "aquamarine2", "deepskyblue3", "gold1", "orange2", "darkorchid")
shapes <- c(15,16,17,18)
scatter_plot<-ggplot(data = mean_scores, aes(x=average,y= R_value, color = as.factor(int), shape = as.factor(outline),label=rownames(mean_scores))) 

# Add the label Milestones
scatter_plot<-scatter_plot+geom_hline(yintercept=0, color = "black", size=0.3)
scatter_plot<-scatter_plot+geom_vline(xintercept=log(1), color = "black", size=0.3)
scatter_plot<-scatter_plot+scale_color_manual("Contacts",values=pointcols, labels=c("I1", "I2_A", "I2_P", "I3_A", "I3_P", "I4", "I5", "I6"), drop = FALSE) 
scatter_plot<-scatter_plot+scale_shape_manual("Scores",values=shapes, labels=c("PCC2", "PCC1", "APCC1", "APCC2"), drop = FALSE) 
scatter_plot<-scatter_plot+geom_point(data = mean_scores , aes(y=R_value, color=as.factor(int), shape = as.factor(outline)), size=4) 
scatter_plot<-scatter_plot+geom_point(data = mean_groups, aes(y = R_value, x= average, shape = as.factor(outline), color = as.factor(int), label = rownames(mean_groups)), size=5) 

#theme the plot
scatter_plot<-scatter_plot+ theme_bw() + 
  ylim(-1,1) +
  theme(panel.background = element_rect(fill = "grey95"),
        text = element_text(family = "sans", colour = "#3E4648", size = 10), 
        plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
        aspect.ratio = 1/1.2,
        plot.margin = unit(c(0,0,0,0), "cm"))

# Print plot
scatter_plot
ggsave(filename = paste("correlation_avarage_", score_name,".svg", sep = ""), path = "E:/DENV_dimer_project/Output/bZIP tests")
