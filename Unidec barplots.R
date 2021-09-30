library(ggplot2)
library(readxl)

Peak_data <- read_excel("E:/native/peptides/Unidec_spectra/Peak_data.xlsx", sheet = "All")
Peak_data <- read_excel("E:/native/peptides/Unidec_spectra/peak_data2.xlsx")
Peak_data$Peptide <- factor(Peak_data$Peptide, levels = c("CC485scr", "CC78", "CC116", "CC233", "CC350", "CC396", "CC485"))
Peak_data$Oligomer <- factor(Peak_data$Oligomer, levels = c("Monomer", "Dimer", "Trimer", "Tetramer", "Pentamer"))


Native_MS_plot <- ggplot(data=Peak_data, aes(x=Peptide, y= as.numeric(Intensity)*100, fill = Oligomer)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values = c("#ed1f24", "#4995d0", "#e8e621", "#91479b", "#e99722"))+
  theme_bw() + 
  theme(panel.background = element_rect(fill = "grey95"),
        text = element_text(family = "sans", colour = "#3E4648", size = 10), 
        plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
        aspect.ratio = 1/1.2,
        plot.margin = grid::unit(c(0,0,0,0), "cm"))+ 
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
  ggtitle("Oligomer intesity of deconvoluted spectra (Unidec)") + 
  labs(colour = "Peptides")
Native_MS_plot 

ggsave("Unidec_intensity_plot2.svg", plot = Native_MS_plot, path = "E:native/peptides")
