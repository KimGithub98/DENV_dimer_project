library(ggplot2)

#convert CD mdeg to MRME
blanc <- read.delim("E:/CD/Kim (Joost)/blanc2.txt", header = FALSE)
colnames(blanc) <- c("Wavelength", "CD", "HT", "Abs")
CC_output_all <- data.frame()
stats_all <- data.frame()

for (CC in c("CC78", "CC116", "CC233", "CC350", "CC396", "CC485", "CC485scr")){
  CC_output <- read.delim(paste("E:/CD/Kim (Joost)/", CC, "_A.txt", sep = ""), header = FALSE)
  colnames(CC_output) <- c("Wavelength", "CD", "HT", "Abs")
  Info = Peptide_inf[Peptide_inf$peptide == CC,]
  MW = Info[3]
  N = Info[2]
  MWR = MW/(N-1)
  c = Info[4]
  MRME_conversion <- function(Oobs, MWR_in, c_in){
    mean_residue_molar_ellipticity = (Oobs*MWR_in)/(10*0.1*c_in)
    return(mean_residue_molar_ellipticity)
  }
  MRME_list = c()
  for (WL in CC_output$Wavelength){
    blanc_sub = blanc$CD[blanc$Wavelength == WL]
    Obs = CC_output$CD[CC_output$Wavelength == WL] - blanc_sub
    MRME <- MRME_conversion(Oobs = Obs, MWR_in = MWR, c_in = c)
    MRME_list = c(MRME_list, MRME$MW)
  }
  group_list = rep(CC, nrow(CC_output))
  CC_output <- cbind(CC_output, MRME_list, group_list)
  colnames(CC_output) <- c("Wavelength", "CD", "HT", "Abs", "MRME", "group")
  
  ##find mimima 208 (200-215) and 222 (210-230)
  min208 = min(CC_output$MRME[CC_output$Wavelength > 200 & CC_output$Wavelength < 215])
  min222 = min(CC_output$MRME[CC_output$Wavelength > 210 & CC_output$Wavelength < 230])
  max214 = max(CC_output$MRME[CC_output$Wavelength > 208 & CC_output$Wavelength < 222])
  ratio = min222/min208
  min222opt = -28500
  helicity = (min222/min222opt)*100
  print(ratio)
  print(helicity)
  stats_all = rbind(stats_all, c(CC, ratio, helicity, as.numeric(max214)))
  CC_output_all = rbind(CC_output_all, CC_output)
}

colnames(stats_all) <- c("CC", "ratio", "helicity", "max214")
stats_all = stats_all[order(stats_all$CC),]

CD_plot <- ggplot(data=CC_output_all, aes(x=Wavelength, y =MRME, group = group)) + 
  geom_line(aes(color = group), size = 1) + 
  scale_color_manual(values = c("chartreuse3", "red4", "brown2", "aquamarine2", "deepskyblue3", "gold1", "orange2", "darkorchid"))+
  annotate(geom = "text", x = 215, y = as.numeric(stats_all$max214), label = paste("ratio:", round(as.numeric(stats_all$ratio), digits = 2), "helicity:", round(as.numeric(stats_all$helicity), digits = 1)), color = c("chartreuse3", "red4", "brown2", "aquamarine2", "deepskyblue3", "gold1", "orange2")) + 
  scale_x_continuous(limits = c(190, 250), expand = c(0,0)) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "grey95"),
        text = element_text(family = "sans", colour = "#3E4648", size = 10), 
        plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
        aspect.ratio = 1/1.2,
        plot.margin = grid::unit(c(0,0,0,0), "cm"))+ 
  ggtitle("Circular Dichroism of top 6 peptides") + 
  labs(colour = "Peptides")

ggsave("CD_plot.svg", plot = CD_plot, path = "E:CD")
