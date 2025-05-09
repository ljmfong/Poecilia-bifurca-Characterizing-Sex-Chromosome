######################################
###### ASE density and boxplot #######
######################################

rm(list=ls())
ls() 

library(ggplot2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(patchwork)
require(grid)
require(cowplot) 

#############################################
####### Make the density plot first #########
#############################################


getwd()
#setwd()
#Set to the appropriate directory

data_f <- read.csv("ASE_GATK/test_mafs/female_medians.csv",stringsAsFactors=F,header=T,sep=",")
data_m <- read.csv("ASE_GATK/iulia_scripts/consistent_ase_male.csv",stringsAsFactors=F,header=T,sep=",")

data_f_auto <- subset(data_f, chromo == "autosome")
data_f_sc <- subset(data_f, chromo == "sexchromo")
data_m_auto <- subset(data_m, chromo == "autosome")
data_m_sc <- subset(data_m, chromo == "sexchromo")

median(data_f_auto$average_majallele_fraction)
#[1] 0.6738552
median(data_f_sc$average_majallele_fraction)
#[1] 0.6917556

median(data_m_auto$average_majallele_fraction)
#[1] 0.6666667
median(data_m_sc$average_majallele_fraction)
#[1] 0.7673183


#Stat Test:


test_data_m_table <- table(data_m$chromo, data_m$average_majallele_fraction)

chisq_test(test_data_m_table)
# n statistic     p    df method          p.signif
#2251     1391. 0.0276  1292 Chi-square test *       
chisq_test(table(data_f$chromo, data_f$average_majallele_fraction))
# n statistic     p    df method          p.signif
#32224    10099.     1 11039 Chi-square test ns   


####### PLOT ########

plot_f = ggplot(data_f, aes(x=average_majallele_fraction, fill=chromo)) + 
  geom_density(alpha=0.7, bw = 0.015, adjust = 1.8) + 
  geom_vline(aes(xintercept=0.6738552),color="slategray3",linetype="dashed",size=0.8) + 
  geom_vline(aes(xintercept=0.6917556),color="seagreen4",linetype="dashed",size=0.8) +
  scale_fill_manual(values=c("slategray3", "seagreen4"), labels = c("Autosomes", "Sex chromosomes")) +
  coord_cartesian(xlim=c(0.5,1)) +
  coord_cartesian(ylim=c(0,5.1)) +
  theme(
    text=element_text(size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    plot.title =element_text(color="black", size=0),
    axis.line.y = element_line(color="black", size = 1),
    axis.line.x = element_line(color="black", size = 1),
    axis.text.x = element_text(size=20,color="black"),
    axis.text.y = element_text(size=20,color="black"),
    axis.title.x=element_text(size=24, color="black"), 
    axis.title.y=element_text(size=24, color="black"),
    legend.position="top",
    legend.text = element_text(colour="black",size=16),
    legend.title= element_text(colour="black",size=0)
  ) +
  scale_x_continuous(limit=c(0.5,1)) +
  xlab('Major Allele Ratio') +
  ylab(expression('Density'))
#annotate("text", x = 0.74, y = 6,5, label = "p = 0.13" , color="black", size=4 , angle=0, fontface="bold")
#geom_text(x=0.74, y=6.5, label="p = 0.13", size=4)

plot_f

plot_m = ggplot(data_m, aes(x=average_majallele_fraction, fill=chromo)) + 
  geom_density(alpha=0.7, bw = 0.015, adjust = 1.8) + 
  geom_vline(aes(xintercept=0.6666667),color="slategray3",linetype="dashed",size=0.8) + 
  geom_vline(aes(xintercept=0.7673183),color="seagreen4",linetype="dashed",size=0.8) +
  scale_fill_manual(values=c("slategray3", "seagreen4"), labels = c("Autosomes", "Sex chromosomes")) +
  coord_cartesian(xlim=c(0.5,1)) +
  coord_cartesian(ylim=c(0,5.1)) +
  theme(
    text=element_text(size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    plot.title =element_text(color="black", size=0),
    axis.line.y = element_line(color="black", size = 1),
    axis.line.x = element_line(color="black", size = 1),
    axis.text.x = element_text(size=20,color="black"),
    axis.text.y = element_text(size=20,color="black"),
    axis.title.x=element_text(size=24, color="black"), 
    axis.title.y=element_text(size=24, color="black"),
    legend.position="top",
    legend.text = element_text(colour="black",size=16),
    legend.title= element_text(colour="black",size=0)
  ) +
  scale_x_continuous(limit=c(0.5,1)) +
  xlab('Major Allele Ratio') +
  ylab(expression('Density')) +
  geom_text(x=0.72, y=5.2, label="*", size=8) + geom_segment(aes(x=0.6666667,y=5.1,xend=0.7673183,yend=5.1), size = 1.3)

plot_m

combo <- plot_f + plot_m + plot_layout(guides = "collect") & theme(legend.position = "top")
combo


#############################################
########## Now make the boxplots ############
#############################################


data = read.table(file.choose(),stringsAsFactors=F,header=T,sep="\t")
#Using files_for_fig1/bifurca_gene_pos_all_rpkm_CLEAN.txt


data_autosomes <- subset(data, LG != "LG12")
data_sexchromo <- subset(data, LG == "LG12")
data_sexchromo$inv_start <- 34233588
data_sexchromo$correct_start <- data_sexchromo$inv_start - data_sexchromo$pos 
data_nonrecombo <- subset(data_sexchromo, correct_start < 30000000)
data_PAR <- subset(data_sexchromo, correct_start > 30000000)
data_autoC <- rbind(data_PAR[,1:12], data_autosomes)

data_auto_male_expr = log2(data_autoC$avg_mal)
data_auto_female_expr = log2(data_autoC$avg_fem)

data_x_male_expr = log2(data_nonrecombo$avg_mal)
data_x_female_expr = log2(data_nonrecombo$avg_fem)



data_autoC$chromo <- "autosome"
data_nonrecombo$chromo <- "nonrecombo"

commcols <- intersect(names(data_autoC), names(data_nonrecombo))
all_data_combo <- bind_rows(
  select(data_autoC, all_of(commcols)),
  select(data_nonrecombo, all_of(commcols))
)

chisq_test(table(all_data_combo$chromo, all_data_combo$avg_mal))
#Output:
n statistic     p    df method          p.signif
* <int>     <dbl> <dbl> <int> <chr>           <chr>   
  1 16639    16554. 0.428 16522 Chi-square test ns  

chisq_test(table(all_data_combo$chromo, all_data_combo$avg_fem))
#Output:
n statistic     p    df method          p.signif
* <int>     <dbl> <dbl> <int> <chr>           <chr>   
  1 16639    16537. 0.303 16444 Chi-square test ns 

par(mar=c(5,6,4,1)+.1)
boxplot(data_auto_female_expr, data_x_female_expr, data_auto_male_expr, data_x_male_expr,
        notch=T, outline=F, ylab=expression('Log'[2]*' Expression'),
        names = c("Female A", "Female X", "Male A", "Male X"),
        border=c("red","red", "blue", "blue"),
        col = NULL,
        boxwex = 0.5, boxlwd = 2.5, whisklwd = 2.5,
        cex.lab=1.8, cex.axis = 1.8, cex.main = 3)

