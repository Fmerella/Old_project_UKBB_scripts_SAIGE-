
library(tidyverse)
library(data.table)

library(ggrepel)
library(readr)

resultschrom7part1<-fread("concatenated_chrom7part1MAC05_results.txt", sep = "\t")

resultschrom7part1$p.value<-as.character(resultschrom7part1$p.value)

class(resultschrom7part1$p.value)

sum(is.na(resultschrom7part1$POS)) ## 0

sum(is.na(resultschrom7part1$BETA))  ## 0
###converts p-values to -log10 p-values(important in the case of extremely small p-vals that may be rounded by R)

resultschrom7part1$logP = -log10(
  as.numeric(
    gsub("E-*","",resultschrom7part1$p.value)
  ) * 0.1
) +
  as.numeric(
    gsub(".*E-", "",resultschrom7part1$p.value)
  ) - 1

head(resultschrom7part1)

sum(is.na(resultschrom7part1$POS))
sum(is.na(resultschrom7part1$logP))### 13
sum(is.na(resultschrom7part1$CHR)) ### 0

na_logP_results<-subset(resultschrom7part1, is.na(logP))

### There are 13 

library(magrittr)
library(ggplot2)

####for point shape to indicate effect direction

man_dfchrom7part1= resultschrom7part1%>%
  arrange(POS) %>%
  mutate(shape = if_else(BETA > 0, "+", "-"),
         shape = factor(shape, levels = c("+","-")))

man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38180264] = "rs7349990"
man_dfchrom7part1$SNP[man_dfchrom7part1$PO  == 38187025] = "rs1784487874"
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38221378] = "rs34938473" ### until here is STARDNL
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38246730] = "rs35616688" #### not inside exon but inside the TRG locus
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38249215] = "rs35633605" ### not inside exon but inside the TRG locus
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38254708] = "rs61181871" #### not inside exon but inside TRG locus
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38261835] = "rs13230963" ### inside the TARP gene
man_dfchrom7part1$SNP[man_dfchrom7part1$POS == 38279290] = "rs35084352"  ### until here is TARP


write_csv(man_dfchrom7part1,"man_dfchrom7part1.csv")


man_dfchrom7part1$shape<-as.character(man_dfchrom7part1$shape)

### Select further based on logP TO find variants that needs to be named on the plot ###


man_dfchrom7part1selected<-man_dfchrom7part1 %>%
  filter(man_dfchrom7part1$logP >5.745)

dim(man_dfchrom7part1selected)  ### 8 26

man_dfchrom7part1selected<-man_dfchrom7part1selected%>%
                    arrange(logP)


man_dfchrom7part1selected<-man_dfchrom7part1selected%>%
  select(CHR,POS,MarkerID,BETA,p.value,logP,shape)

#this can be used to annotate genes, please comment out as necessary
## gene1 = c(38178245-38230670)
## gene2 = c(38240024-38273636)


write_csv(man_dfchrom7part1selected,"man_dfchrom7part1selected.csv")



sum(is.na(man_dfchrom7part1$BETA))  ## 14

man_dfchrom7part1$BETA<-as.numeric(man_dfchrom7part1$BETA)
man_dfchrom7part1$POS<-as.numeric(man_dfchrom7part1$POS)


graph_title = 'Manhattan Plot for the Region 38079029-38339030 of Chromosome 7 part1'


file_name = "Saigechrom7part1_results.png"


options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)


ggplot(man_dfchrom7part1, aes(x = POS, y = logP)) +
  geom_point(aes(shape = shape, fill = BETA), colour = "grey50", size = 2) +
  scale_shape_manual(values = c(25, 24), name = "Beta Sign") +
  scale_fill_gradient2(low = "blue", mid = "lightyellow1", midpoint = 0, high = "red", name = "Beta") +
  
  ylim(-3, 10)+
  
  # Repelling SNP labels
  geom_text_repel(aes(label = SNP), size = 3, box.padding = 0.5,
                  data = man_dfchrom7part1[man_dfchrom7part1$SNP %in% c("rs7349990", "rs1784487874", "rs34938473", "rs35616688", "rs35633605", "rs61181871", "rs13230963", "rs35084352"),]) +
  
  theme(axis.text.x = element_blank()) +
  xlab("Position") +
  ylab("log10(p)") +
  
  ### Drawing the green rectangles ###
  
  annotate("rect", xmin = 38178245, xmax = 38230670, ymin = -0.005, ymax = 0, alpha= 0.75, fill = "darkgreen") +
  annotate("rect", xmin = 38240024, xmax = 38368055, ymin = -0.005, ymax = 0, alpha= 0.75, fill = "darkgreen") +
  annotate("text", x = 38195000, y = -0.01, label = "STARD3NL", size = 4, color ="white", fontface ="bold", hjust = 0.5) +
  annotate("text", x = 38259643, y = -0.01, label = "TARP", size = 4, color ="yellow", fontface = "bold") +
  
  geom_hline(yintercept = 7.30103, linetype = "dotted", color = 'red', size = 0.5) +
  
  ggtitle(graph_title) +
  ylab('log10 p-value') +
  
  
  xlim(min(man_dfchrom7part1$POS) - 5.0 * diff(range(man_dfchrom7part1$POS)),
       max(man_dfchrom7part1$POS) + 5.0 * diff(range(man_dfchrom7part1$POS)))



ggsave(file_name, dpi = 300, height = 7, width = 20)


### Alternative codes using Cairo

install.packages("Cairo")
library(Cairo)

# Load required libraries
library(ggplot2)
library(ggrepel)

# Save plot as a PNG using Cairo, which does not need X11

CairoPNG("ManPlot_chrom7part1.png", width = 800, height = 600)
  
  # Create the Manhattan plot
  ggplot(man_dfchrom7part1, aes(x = POS, y = logP)) +
    geom_point(aes(shape = shape, fill = BETA), colour = "grey50", size = 2) +
    scale_shape_manual(values = c(25, 24), name = "Beta Sign") +
    scale_fill_gradient2(low = "blue", mid = "lightyellow1", midpoint = 0, high = "red", name = "Beta") +
    ylim(-3, 10) +
    
    geom_text_repel(aes(label = SNP), size = 3, box.padding = 0.5,
                    data = man_dfchrom7part1[man_dfchrom7part1$SNP %in% c("rs7349990", "rs1784487874", "rs34938473", "rs35616688", "rs35633605", "rs61181871", "rs13230963", "rs35084352"),]) +
    
    theme(axis.text.x = element_blank()) +
    xlab("Position") +
    ylab("log10(p)") +
    
    annotate("rect", xmin = 38178245, xmax = 38230670, ymin = -0.005, ymax = 0, alpha= 0.75, fill = "darkgreen") +
    annotate("rect", xmin = 38240024, xmax = 38368055, ymin = -0.005, ymax = 0, alpha= 0.75, fill = "darkgreen") +
    annotate("text", x = 38195000, y = -0.01, label = "STARD3NL", size = 4, color = "white", fontface = "bold", hjust = 0.5) +
    annotate("text", x = 38259643, y = -0.01, label = "TARP", size = 4, color = "yellow", fontface = "bold") +
    
    geom_hline(yintercept = 7.30103, linetype = "dotted", color = 'red', size = 0.5) +
    ggtitle(graph_title) +
    ylab('log10 p-value') +
    xlim(min(man_dfchrom7part1$POS) - 5.0 * diff(range(man_dfchrom7part1$POS)),
         max(man_dfchrom7part1$POS) + 5.0 * diff(range(man_dfchrom7part1$POS)))
  
  # Close the PNG device
  dev.off()














