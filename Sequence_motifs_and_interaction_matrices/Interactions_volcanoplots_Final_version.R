# Script title: Interactions_volcanoplots
# Version title: Final_version
# Short application: compares interaction probability matrices and visualizes these in volcanoplots
# Input files: Uses interaction probability and count matrices from dataset and the expected interaction probability matrices
# Output files: .png and .svg file of the volcanoplots in Interaction_volcanoplots and .csv file with the relative probabilities and P_values plotted in these volcanoplots
# How to use: run script after collecting all necessary interaction probability files

# Name: Kim N. Wijnant
# Group: Biomolecular Mass Spectrometry and Proteomics
# Date of last edit: 27/09/2021

library(ggrepel)
library(ggplot2)

##plots interaction frequency data as volcano plots

P_values <- data.frame()

## statistical testing with multinomial distribution
stat_test <- function(input_1, input_2_P){
  files_input_1 <- read.csv(input_1)
  files_input_2 <- read.csv(input_2_P)
  rownames(files_input_1) <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")
  rownames(files_input_2) <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")
  for(AA_in in colnames(files_input_1)[1:(length(colnames(files_input_1))-1)]){
    for(int_in in rownames(files_input_1)){
      sum_1 <- files_input_1[int_in, AA_in]
      sum_int_1 <- sum(files_input_1[int_in, 1:(length(colnames(files_input_1))-1)])
      sum_2 <- files_input_2[int_in, AA_in]
      sum_int_2 <- sum(files_input_2[int_in,])
      probability = sum_2
      if(probability == 0){
        probability = min(files_input_2[int_in,][files_input_2[int_in,] != 0])
      }
      tt = dmultinom(c(sum_1, (sum_int_1-sum_1)), prob = c(probability, 1-probability))
      P_values[int_in, AA_in] <- tt
    }
  }
  return(P_values)
}

##function for statistical testing multinomial distribution with expected file
stat_test_expected <- function(input_1, input_2_P){
  files_input_1 <- read.csv(input_1)
  files_input_2 <- input_2_P
  rownames(files_input_1) <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")
  rownames(files_input_2) <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")
  for(int_in in rownames(files_input_1)){
    sum_int_1 <- sum(files_input_1[int_in, 1:(length(colnames(files_input_1))-1)])
    new_files_input_2 <- round(files_input_2[int_in,]*sum_int_1)
    sum_int_2 <- sum(new_files_input_2)
    for(AA_in in colnames(files_input_1)[1:(length(colnames(files_input_1))-1)]){
      sum_1 <- files_input_1[int_in, AA_in]
      sum_int_1 <- sum(files_input_1[int_in, 1:(length(colnames(files_input_1))-1)])
      sum_2 <- new_files_input_2[int_in, AA_in]
      probability = sum_2/sum_int_2
      if(probability == 0){
        probability = 1/sum_int_2
      }
      tt = dmultinom(c(sum_1, (sum_int_1-sum_1)), prob = c(probability, 1-probability))
      P_values[int_in, AA_in] <- tt
    }
  }
  return(P_values)
}

#get name of a variable
get_var_name <- function(v1) {
  deparse(substitute(v1))
}

#function for PCC volcano plot
table_volcano_plot_PCC <- function(Probibility_table, P_value_table){
  data2 <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(data2) <- c("AA", "Int", "Probibility", "P_value")
  for(col in colnames(Probibility_table)){
    for(row in c(1:4,5:8)){
      Prob <- Probibility_table[row, col]
      P_value <- P_value_table[row, col]
      row_name <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")[row] #combines I4/1 as these are the same in PCC
      new_row <- c(col, row_name, as.numeric(Prob), as.numeric(P_value))
      data2 <- rbind(data2, new_row)
    }
  }
  colnames(data2) <- c("AA", "Int", "Probibility", "P_value")
  return(data2)
}

#function for APCC volcano plot
table_volcano_plot_APCC <- function(Probibility_table, P_value_table){
  print(Probibility_table[,1:4])
  data2 <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(data2) <- c("AA", "Int", "Probibility", "P_value")
  for(col in colnames(Probibility_table)){
    for(row in c(1:4, 5:8)){
      Prob <- Probibility_table[row, col]
      P_value <- P_value_table[row, col]
      row_name <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")[row] #combines I2/3 as these interactions are the same in APCC
      new_row <- c(col, row_name, as.numeric(Prob), as.numeric(P_value))
      data2 <- rbind(data2, new_row)
    }
  }
  colnames(data2) <- c("AA", "Int", "Probibility", "P_value")
  return(data2)
}

#creates csv and volcano plots for 
for (database in c("CCplus_")){
  i_csv = 0
  for (orientation in c("APCC", "PCC")){ #makes plots for APCC and PCC projections
    #get expected interaction matrices for the orientation
    if (orientation == "APCC"){
      interaction_frequency_CC = read.csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/Expected_interaction_matrix_APCC_24_28.csv", sep=",", row.names = 1)
      point_col = c("chartreuse3", "brown2", NA, "deepskyblue3", NA, "gold1", "orange2", "darkorchid") #colors APCC
    }
    else if (orientation == "PCC"){
      interaction_frequency_CC = read.csv("E:/DENV_dimer_project/Output/Expected_interaction_matrix/Expected_interaction_matrix_PCC_24_28.csv", sep=",", row.names = 1)
      point_col = c("chartreuse3", NA, "brown2", NA,"deepskyblue3", "gold1", "orange2", "darkorchid") #colors PCC
    }
    # get interaction probability matrices for orientation projected as orientation and other orientation projected as orientation
    file1 <- paste("E:/DENV_dimer_project/Output/Make_interaction_matrices/Count/", database, orientation, "_length_24_28_cdhit_filtered/as_", orientation, ".csv", sep = "")
    if (orientation == "APCC"){other_orientation = "PCC"}
    else if (orientation == "PCC"){other_orientation = "APCC"}
    file2 <- paste("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/", database, other_orientation, "_length_24_28_cdhit_filtered/as_", orientation, "_P.csv", sep = "")
    input_1 <- read.csv(paste("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/", database, orientation, "_length_24_28_cdhit_filtered/as_", orientation, "_P.csv", sep = ""))
    input_2 <- read.csv(file2)
    input_3 <- interaction_frequency_CC
    input_3$Sum <- c(1,1,1,1,1,1,1,1) #adds sum column to make dataframes equal length
    
    #calculate p values
    P_values_APCCvsPCC <- stat_test(file1, file2)
    P_values_expected <- stat_test_expected(file1, interaction_frequency_CC)
    
    #calculate relative probabilities
    rel_prob_APCCvsPCC <- input_1-input_2
    rel_prob_expected <- input_1-input_3
    
    #get data for volcanoplots
    if (orientation == "APCC"){
      data_APCCvsPCC <- table_volcano_plot_APCC(rel_prob_APCCvsPCC, P_values_APCCvsPCC)
      data_expected <- table_volcano_plot_APCC(rel_prob_expected, P_values_expected)
    }
    if (orientation == "PCC"){
      data_APCCvsPCC <- table_volcano_plot_PCC(rel_prob_APCCvsPCC, P_values_APCCvsPCC)
      data_expected <- table_volcano_plot_PCC(rel_prob_expected, P_values_expected)
    }
    
    #only visualize labels of P>0.5 interactions
    data_APCCvsPCC$delabel[as.numeric(data_APCCvsPCC$P_value) < 0.05 & as.numeric(data_APCCvsPCC$Probibility) < -0.02] <- data_APCCvsPCC$AA[as.numeric(data_APCCvsPCC$P_value) < 0.05 & as.numeric(data_APCCvsPCC$Probibility) < -0.02] 
    data_APCCvsPCC$delabel[as.numeric(data_APCCvsPCC$P_value) < 0.05 & as.numeric(data_APCCvsPCC$Probibility) > 0.02] <- data_APCCvsPCC$AA[as.numeric(data_APCCvsPCC$P_value) < 0.05 & as.numeric(data_APCCvsPCC$Probibility) > 0.02]
    data_expected$delabel[as.numeric(data_expected$P_value) < 0.05 & as.numeric(data_expected$Probibility) < -0.02] <- data_expected$AA[as.numeric(data_expected$P_value) < 0.05 & as.numeric(data_expected$Probibility) < -0.02] 
    data_expected$delabel[as.numeric(data_expected$P_value) < 0.05 & as.numeric(data_expected$Probibility) > 0.02] <- data_expected$AA[as.numeric(data_expected$P_value) < 0.05 & as.numeric(data_expected$Probibility) > 0.02]
    
    #create volcano plots
    i= 0
    for(data in list(data_APCCvsPCC, data_expected)){
      if(i == 0){data_name = "APCCvsPCC_"}
      if (i == 1){data_name = "expected_"}
      
      volcano <- ggplot(data = data, aes(x = as.numeric(Probibility), y = -log10(as.numeric(P_value)), color = as.factor(Int), label = delabel)) +
        geom_point() +
        scale_color_manual(values= point_col) + 
        geom_vline(xintercept=c(-0.02, 0.02), col="grey70", lty= "dashed") + 
        geom_hline(yintercept=-log10(0.05), col="grey70", lty= "dashed") + 
        geom_text_repel() + 
        xlab("relative probability") + 
        ylab(expression("-Log"[10]*"P")) + 
        theme_bw() + 
        theme(panel.background = element_rect(fill = "grey95"),
              text = element_text(family = "sans", colour = "#3E4648", size = 10), 
              plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
              aspect.ratio = 1/1.2,
              plot.margin = unit(c(0,0,0,0), "cm"))+ 
        ggtitle(paste(database, data_name, "as", orientation, sep = "")) + 
        labs(colour = "CC interactions")
      
      #save volcanoplots
      ggsave(paste(database, data_name, "as", orientation, ".png", sep = ""), plot = volcano, path = "E:/DENV_dimer_project/Output/Interaction_volcanoplots")
      ggsave(paste(database, data_name, "as", orientation, ".svg", sep = ""), plot = volcano, path = "E:/DENV_dimer_project/Output/Interaction_volcanoplots")
      write.csv(data, paste("E:/DENV_dimer_project/Output/Interaction_volcanoplots/", database, data_name, "as", orientation, ".csv", sep = ""), row.names = FALSE)
      
      #saves relative probabilities in a dataframe
      if (i_csv == 0){
        R_interactions = cbind(data$AA, data$Int, data$Probibility)
      }
      else if (i_csv >= 0){
        R_interactions = cbind(R_interactions, data$Probibility)
      }
      i_csv = i_csv + 1
      i = i + 1
    }
    }
}

#save Relative probability 
colnames(R_interactions) <- c("AA", "Int", "CCplus_APCCvsPCC_asAPCC", "CCplus_expected_asAPCC", "CCplus_APCCvsPCC_asPCC", "CCplus_expected_asPCC")
write.csv(R_interactions, "E:/DENV_dimer_project/Output/Interaction_volcanoplots/rel_prob_for_every_aa_int.csv")

#save non-relative probability
table_non_relative_PCC <- function(Probibility_table){
  data2 <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(data2) <- c("AA", "Int", "Probibility")
  for(col in colnames(Probibility_table)){
    for(row in 1:8){
      Prob <- Probibility_table[row, col]
      row_name <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")[row] #combines I4/1 as these are the same in PCC
      new_row <- c(col, as.character(row_name), as.numeric(Prob))
      data2 <- rbind(data2, new_row)
    }
  }
  colnames(data2) <- c("AA", "Int", "Probibility")
  return(data2)
}
table_non_relative_APCC <- function(Probibility_table){
  data2 <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(data2) <- c("AA", "Int", "Probibility")
  for(col in colnames(Probibility_table)){
    for(row in 1:8){
      Prob <- Probibility_table[row, col]
      row_name <- c("I2_P", "I3_P", "I2_A", "I3_A", "I1", "I4", "I5", "I6")[row] #combines I2/3 as these interactions are the same in APCC
      new_row <- c(col, as.character(row_name), as.numeric(Prob))
      #print(new_row)
      data2 <- rbind(data2, new_row)
    }
  }
  colnames(data2) <- c("AA", "Int", "Probibility")
  return(data2)
}

CCplus_APCC_asAPCC <- read.csv("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/CCplus_APCC_length_24_28_cdhit_filtered/as_APCC_P.csv", na.strings = "NNN")
CCplus_PCC_asAPCC <- read.csv("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/CCplus_PCC_length_24_28_cdhit_filtered/as_APCC_P.csv", na.strings = "NNN")
CCplus_APCC_asPCC <- read.csv("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/CCplus_APCC_length_24_28_cdhit_filtered/as_PCC_P.csv", na.strings = "NNN")
CCplus_PCC_asPCC <- read.csv("E:/DENV_dimer_project/Output/Make_interaction_matrices/Probability/CCplus_PCC_length_24_28_cdhit_filtered/as_PCC_P.csv", na.strings = "NNN")

CCplus_APCC_asAPCC <-table_non_relative_APCC(CCplus_APCC_asAPCC)
CCplus_PCC_asAPCC <-table_non_relative_APCC(CCplus_PCC_asAPCC)
CCplus_APCC_asPCC <-table_non_relative_PCC(CCplus_APCC_asPCC)
CCplus_PCC_asPCC <-table_non_relative_PCC(CCplus_PCC_asPCC)

CCplus <- cbind(CCplus_APCC_asAPCC$AA, CCplus_APCC_asAPCC$Int, CCplus_APCC_asAPCC$Probibility, CCplus_PCC_asAPCC$Probibility, CCplus_APCC_asPCC$Probibility, CCplus_PCC_asPCC$Probibility)
colnames(CCplus) <- c("AA", "Int", "CCplus_APCC_asAPCC", "CCplus_PCC_asAPCC", "CCplus_APCC_asPCC", "CCplus_PCC_asPCC")
write.csv(CCplus, "E:/DENV_dimer_project/Output/Interaction_volcanoplots/NR_prob_for_every_aa_int.csv")

