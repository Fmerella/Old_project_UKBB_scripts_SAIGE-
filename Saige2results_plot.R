library(tidyverse)
library(data.table)
library(ggrepel)
library(readr)
getwd() ## "/home/rstudio-server/Saige2chrom7"


#read in concatenated results files

resultschrom7<-fread("concatenated_results.txt", sep = "\t")

resultschrom7$p.value<-as.character(resultschrom7$p.value)

class(resultschrom7$p.value)


#converts p-values to -log10 p-values (this is particularly important in the case of extremely small p-vals that may be rounded by R)

resultschrom7$logP = -log10(
  as.numeric(
    gsub("E-*","",resultschrom7$p.value)
  ) * 0.1
) +
  as.numeric(
    gsub(".*E-", "",resultschrom7$p.value)
  ) - 1

head(resultschrom7)

library(magrittr)
library(ggplot2)

####for point shape to indicate effect direction
man_dfchrom7 = resultschrom7%>%
  arrange(POS) %>%
  mutate(shape = if_else(BETA > 0, "+", "-"),
         shape = factor(shape, levels = c("+","-")))


graph_title = 'Manhattan Plot for the Region 38079029 - 38339030 of Chromosome 7'

library(ggrepel)

#name of output file, should end with png or pdf dependent on image type wanted
file_name = "Saigechrom7_results.png"

man_dfchrom7<-as.character(man_dfchrom7$shape)

### Select further based on logP TO find variants that needs to be named on the plot ###

man_dfchrom7selected<-man_dfchrom7 %>%
                       filter(man_dfchrom7$logP >5.745)


man_dfchrom7selected$logP


#### #optional - add labels for particular points. Comment out if not using. Can also label with gene names etc
### man_df$SNP[man_df$POS == 38187025] = "rs1583773019" ### connected to gene STARD3NL
### man_df$SNP[man_df$POS == 38221378] = "rs34938473"   ### connected to gene STARD3NL
### man_df$SNP[man_df$POS == 38261835] = "rs13230963"  ### connected to TARP gene
### man_df$SNP[man_df$POS == 38279290] = "rs35084352"



man_dfchrom7$SNP[man_dfchrom7$POS == 38187025] = "rs1583773019"
man_dfchrom7$SNP[man_dfchrom7$POS == 38221378] = "rs34938473"   
man_dfchrom7$SNP[man_dfchrom7$POS == 38261835] = "rs13230963"  
man_dfchrom7$SNP[man_dfchrom7$POS == 38279290] = "rs35084352"


#this can be used to annotate genes, please comment out as necessary
### gene1 = c(38178245-38230670)
### gene2 = c(38259643-38273636)


library(ggrepel)



#read in concatenated results files

resultschrom7<-fread("concatenated_results.txt", sep = "\t")

resultschrom7$p.value<-as.character(resultschrom7$p.value)

class(resultschrom7$p.value)


#converts p-values to -log10 p-values (this is particularly important in the case of extremely small p-vals that may be rounded by R)

resultschrom7$logP = -log10(
  as.numeric(
    gsub("E-*","",resultschrom7$p.value)
  ) * 0.1
) +
  as.numeric(
    gsub(".*E-", "",resultschrom7$p.value)
  ) - 1

head(resultschrom7)

library(magrittr)
library(ggplot2)


#name of output file, should end with png or pdf dependent on image type wanted
file_name = "Saigechrom7_results.png"

man_dfchrom7<-as.character(man_dfchrom7$shape)

### Select further based on logP TO find variants that needs to be named on the plot ###

man_dfchrom7selected<-man_dfchrom7 %>%
                       filter(man_dfchrom7$logP >5.745)


man_dfchrom7selected$logP


#### #optional - add labels for particular points. Comment out if not using. Can also label with gene names etc
### man_df$SNP[man_df$POS == 38187025] = "rs1583773019" ### connected to gene STARD3NL
### man_df$SNP[man_df$POS == 38221378] = "rs34938473"   ### connected to gene STARD3NL
### man_df$SNP[man_df$POS == 38261835] = "rs13230963"  ### connected to TARP gene
### man_df$SNP[man_df$POS == 38279290] = "rs35084352"



man_dfchrom7$SNP[man_dfchrom7$POS == 38187025] = "rs1583773019"
man_dfchrom7$SNP[man_dfchrom7$POS == 38221378] = "rs34938473"   
man_dfchrom7$SNP[man_dfchrom7$POS == 38261835] = "rs13230963"  
man_dfchrom7$SNP[man_dfchrom7$POS == 38279290] = "rs35084352"


#this can be used to annotate genes, please comment out as necessary
gene1 = c(38178245-38230670)  ### STARD3NL
gene2 = c(38259643-38273636)  ### TARP GENE

sum(is.na(man_dfchrom7$BETA)) ### 0

man_dfchrom7$BETA<-as.numeric(man_dfchrom7$BETA)

# Assuming `man_dfchrom7` is your dataframe containing necessary SNP and logP data
library(ggrepel)


man_dfchrom7$POS<-as.numeric(man_dfchrom7$POS)




### another alternative version ###

ggplot(man_dfchrom7, aes(x = POS, y = logP)) +
  geom_point(aes(shape = shape, fill = BETA), colour = "grey50", size = 2) +
  scale_shape_manual(values = c(24, 25), name = "Beta Sign") +
  scale_fill_gradient2(low = "blue", mid = "lightyellow1", midpoint = 0, high = "red", name = "Beta") +
  ylim(0, max(man_dfchrom7$logP) + 0.05 * max(man_dfchrom7$logP)) +
  
  # Repelling SNP labels
  geom_text_repel(aes(label = SNP), size = 3, box.padding = 0.5,
                  data = man_dfchrom7[man_dfchrom7$SNP %in% c("rs1583773019", "rs34938473", "rs13230963", "rs35084352"),]) +
  
  theme(axis.text.x = element_blank()) +
  xlab("Position") + ylab("log10(p)") +
  
  annotate("rect", xmin = 38178245, xmax = 38230670, ymin = -10, ymax = 6, alpha = 0.75, fill = "darkgreen") +
  annotate("rect", xmin = 38259643, xmax = 38273636, ymin = -10, ymax = 6, alpha = 0.75, fill = "darkgreen") +
  
  # Repelling gene names
  geom_text_repel(data = data.frame(
    x = c(38195000, 38265000),
    y = c(0, 0),
    label = c("STARD3NL", "TARP GENE")
  ),
  aes(x = x, y = y, label = label),
  size = 3, box.padding = 0.5, point.padding = 0.5) +
  
  geom_hline(yintercept = 7.30103, linetype = "dotted", color = 'red', size = 1) +
  
  ggtitle(graph_title) +
  ylab('log10 p-value') +
  ylim(-10, max(man_dfchrom7$logP)) +
  xlim(min(man_dfchrom7$POS) - 6.0 * diff(range(man_dfchrom7$POS)),
       max(man_dfchrom7$POS) + 6.0 * diff(range(man_dfchrom7$POS)))

ggsave(file_name, dpi = 300, height = 7, width = 18)

