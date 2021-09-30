library("ggplot2")
library("ggpubr")
library(svglite)

## Insert your location of the excel file and the sheets you want to analyse
Location_raw_excel_file = "E:/native/peptides/raw data peptide native MS.xlsx"
Save_location_figures = "E:/native/peptides/raw_spectra" #choose an existing path (this script can not make a new folder)
excel_sheets = c("cc485", "CC78", "CC350", "CC116", "CC233", "CC396", "CC485scr")
lower_limit = 500 #change to lowest M/Z setting
upper_limit = 3000 #change to highest M/Z setting


plot_list <- vector("list", 7)
i = 0
for (excel_sheet in excel_sheets){
  svg(paste(Save_location_figures, "/", excel_sheet, ".svg", sep = ""))
  i = i + 1
  print(excel_sheet)
  print(i)
  raw_data_native_MS <- read_excel(Location_raw_excel_file, sheet = excel_sheet, col_names = FALSE)
  colnames(raw_data_native_MS) <- c("M_Z", "Intensity")
  p <- ggplot(data = raw_data_native_MS, aes(x= M_Z, y = Intensity)) 
  p <- p + geom_line(size = 0.2)
  p <- p + xlab("") + 
    ylab("") +
    ggtitle(excel_sheet)
  p <- p + theme_classic() + 
    #xlim(-0.05,0.05) +
    #ylim(0,4) +
    theme(#panel.background = element_rect(fill = "grey95"),
          #panel.border = element_rect(colour = "#3E4648", fill = NA, size = 1),
          text = element_text(family = "sans", colour = "#3E4648", size = 10), 
          plot.title = element_text(family = "sans", colour = "#3E4648", size = 14, hjust = 0.5),
          aspect.ratio = 1/2,
          plot.margin = unit(c(0,0,0,0), "cm"))  + 
    scale_x_continuous(limits = c(lower_limit, upper_limit), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0 - max(raw_data_native_MS["Intensity"])*0.02, max(raw_data_native_MS["Intensity"])*1.05), expand = c(0,0))
  print(p)
  dev.off()
}



