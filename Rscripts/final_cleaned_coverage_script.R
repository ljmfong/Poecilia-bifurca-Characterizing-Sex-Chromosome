######################################
####### M:F Coverage Analysis ########
######################################

rm(list=ls())
ls() 

library(ggplot2)
library(gridExtra)
library(RcppRoll)
library(zoo)
library(gridExtra)
library(patchwork)
library(tidyverse)


######################################
############ P. bifurca ##############
######################################


#define a moving average function
movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill = NA)
}

windowsizecov = 40

#get data
datacov <- read.csv(file.choose(), stringsAsFactors = FALSE, header=T) 
#This should also have positional info: File is MF_log_avg_cov_pos_23scaffs.csv


datacov$Mlogaverage <- log(datacov$Mcovavg)
datacov$Flogaverage <- log(datacov$Fcovavg)
datacov$MFLogaverage <- (datacov$Mlogaverage)/(datacov$Flogaverage)
datacov$MFlog2avg <- log2((datacov$Mcovavg)/(datacov$Fcovavg))
colnames(datacov)[1] <- "Scaffold"

#Rename the ragtag scaffold chromosomes to P. picta reference for comparison
datacov <- datacov %>%
  mutate(chromosome = case_when(
    Scaffold == "Scaffold_1_RagTag" ~ "LG2",
    Scaffold == "Scaffold_2_RagTag" ~ "LG3",
    Scaffold == "Scaffold_3_RagTag" ~ "LG9",
    Scaffold == "Scaffold_4_RagTag" ~ "LG5",
    Scaffold == "Scaffold_5_RagTag" ~ "LG16",
    Scaffold == "Scaffold_6_RagTag" ~ "LG7",
    Scaffold == "Scaffold_7_RagTag" ~ "LG13",
    Scaffold == "Scaffold_8_RagTag" ~ "LG1",
    Scaffold == "Scaffold_9_RagTag" ~ "LG10",
    Scaffold == "Scaffold_10_RagTag" ~ "LG4",
    Scaffold == "Scaffold_11_RagTag" ~ "LG17",
    Scaffold == "Scaffold_12_RagTag" ~ "LG6",
    Scaffold == "Scaffold_13_RagTag" ~ "LG15",
    Scaffold == "Scaffold_14_RagTag" ~ "LG12",
    Scaffold == "Scaffold_15_RagTag" ~ "LG14",
    Scaffold == "Scaffold_16_RagTag" ~ "LG11",
    Scaffold == "Scaffold_17_RagTag" ~ "LG18",
    Scaffold == "Scaffold_18_RagTag" ~ "LG22",
    Scaffold == "Scaffold_19_RagTag" ~ "LG20",
    Scaffold == "Scaffold_20_RagTag" ~ "LG8",
    Scaffold == "Scaffold_21_RagTag" ~ "LG21",
    Scaffold == "Scaffold_22_RagTag" ~ "LG19",
    Scaffold == "Scaffold_23_RagTag" ~ "LG23"
  ))


autocov <- datacov %>%
  filter_all(any_vars(datacov$Scaffold != "Scaffold_14_RagTag"))

dim(datacov)
#[1] 16518  10
dim(autocov)
#[1] 15833  10

names(datacov)
#[1] "Scaffold"     "binstart"     "binend"       "Fcovavg"      "Mcovavg"      "Mlogaverage"  "Flogaverage"  "MFlog2avg"   
#[9] "chromosome"   "MFLogaverage"

mean(autocov$MFlog2avg)
#[1] 0.0917997

mean(autocov$MFLogaverage)
#[1]  1.021932


#sort data - in this case, the chromosome is already sorted into the bins of 50kB
autocovsorted <- autocov[order(autocov$Scaffold,autocov$binstart),]
dim(autocovsorted)
# [1] 15833     10

#calculate confidence interval for moving average

MFautopermute <- replicate(1000,mean(sample(datacov$MFlog2avg,windowsizecov,replace = F)))
MFautoI25cov <- quantile(MFautopermute, c(.025, .5, .975))[[1]]
MFautoI25cov
#[1] -0.04139033

MFautoI975cov <- quantile(MFautopermute, c(.025, .5, .975))[[3]]
MFautoI975cov
#[1] 0.1320388

#Calculate confidence interval for MFLogaverage for FST gene duplication:

MFautopermute_FST <- replicate(1000,mean(sample(datacov$MFLogaverage,windowsizecov,replace = F)))
MFautoI25cov_FST <- quantile(MFautopermute_FST, c(.025, .5, .975))[[1]]
MFautoI25cov_FST
#[1] 0.9741286

MFautoI975cov_FST <- quantile(MFautopermute_FST, c(.025, .5, .975))[[3]]
MFautoI975cov_FST
#[1] 1.053978


datacov$chromosome <- factor(datacov$chromosome, paste0("LG",1:23), paste0("LG",1:23))

for (i in levels(datacov$chromosome)){
  chr_tester <- subset(datacov, chromosome == i)
  assign(paste0(i), chr_tester)
}


# List of objects starting with "LG"
data_list <- ls(pattern = "^LG")


#LG12 also has been inverted and needs to be flipped:

data12_edit <- LG12
data12_edit$inv_start <- 34233588
data12_edit$correct_start <- data12_edit$inv_start - data12_edit$binstart
data12_edit$binstart <- data12_edit$correct_start
data12_edit <- data12_edit[1:10]
LG12 <- data12_edit



generate_plot <- function(obj) {
  # Perform the operations
  smoothlinefc = movingaverage(obj$MFlog2avg,windowsizecov)
  smoothlinefem = movingaverage(obj$Flogaverage,windowsizecov)
  smoothlinemal = movingaverage(obj$Mlogaverage,windowsizecov)
  
  paneldata = as.data.frame(smoothlinefc)
  paneldata$startmb = obj$binstart
  paneldata$female_expr <- smoothlinefem
  paneldata$male_expr <- smoothlinemal
  paneldata <- na.omit(paneldata)
  
  # Create the plot using ggplot
  plot <- ggplot(obj, aes(x = binstart/1000000, y = MFlog2avg)) +
    geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey70", alpha=0.08)+
    geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
    geom_line(data = paneldata, aes(x= startmb/1000000, y = smoothlinefc), size=1.5) +
    coord_cartesian(ylim=c(-1.2, 1.2)) +
    xlab(NULL) +
    ylab(NULL) + scale_y_continuous(breaks = c(1.0, 0.0, -1.0)) +
    labs(title = obj$chromosome) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 15), axis.text = element_text(size = 14),
          legend.position = 'top')
  
  return(plot)
}

# Use lapply to generate plots for each object in the list
plot_list <- lapply(data_list, function(obj_name) {
  # Access the object
  obj <- get(obj_name)
  
  # Generate the plot
  generate_plot(obj)
})

plot_list_reordered <- plot_list[c(1:1, 12:12, 17:23, 2:11, 13:16)]  

# Create the multi-panel plot using patchwork

multi_panel_plot <- wrap_plots(plotlist = plot_list_reordered, ncol = 4) +
  guide_area() +
  plot_layout(guides = 'collect')
final_plot <- grid.arrange(patchworkGrob(multi_panel_plot), 
                           left = textGrob(expression('M:F log'[2]*' Coverage'), rot = 90, gp = gpar(fontsize = 17)),
                           bottom = textGrob("Chromosome Length (Mb)", gp = gpar(fontsize = 17)))


final_plot



#Plotting Scaffold_14/LG12

datacovScaffold14 <- datacov %>%
  filter_all(any_vars(datacov$Scaffold == "Scaffold_14_RagTag"))

datacovScaffold14 <- datacovScaffold14[order(datacovScaffold14$binstart,as.numeric(levels(datacovScaffold14$binstart))[datacovScaffold14$binstart]),]

dim(datacovScaffold14)
# [1] 685   10

mean(datacovScaffold14$MFlog2avg)
# [1] -0.7815422


mean(datacovScaffold14$MFLogaverage)
#[1] -0.5417238

datacovScaffold14_edit <- datacovScaffold14
datacovScaffold14_edit$inv_start <- 34233588
datacovScaffold14_edit$correct_start <- datacovScaffold14_edit$inv_start - datacovScaffold14_edit$binstart
datacovScaffold14_edit$binstart <- datacovScaffold14_edit$correct_start
datacovScaffold14_edit <- data12_edit[1:10]
datacovScaffold14 <- datacovScaffold14_edit


smoothlinefccov = movingaverage(datacovScaffold14$MFlog2avg, windowsizecov)
smoothlinefemcov = movingaverage(datacovScaffold14$Flogaverage, windowsizecov)
smoothlinemalcov = movingaverage(datacovScaffold14$Mlogaverage, windowsizecov)


paneldatacov <- as.data.frame(smoothlinefccov)
paneldatacov$startmb <- datacovScaffold14$binstart
paneldatacov$females <- smoothlinefemcov
paneldatacov$males <- smoothlinemalcov
names(paneldatacov)
# "smoothlinefccov" "startmb"         "females"         "males"
paneldatacov <- na.omit(paneldatacov)

# plots (you might have to change the axes t o view your plot better - depends on the values plotted)
require(grid)



plotScaffold14 <- ggplot(datacovScaffold14, aes(x= binstart/1000000, y= MFlog2avg)) +
  annotate("rect", xmax = 35, xmin = 30, ymax = 1.2, ymin = -1.2, fill = "forestgreen", alpha = 0.3) +
  geom_rect(xmax=35,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
  geom_line(data = paneldatacov, aes(x=  paneldatacov$startmb/1000000, y= paneldatacov$smoothlinefccov), colour = "black", size=1.5) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,35)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +
  ylim(-1.2,1.2) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Coverage'))

plotScaffold14


######################################
############# P. picta ###############
######################################


datacov_pic <- read.table(file.choose(), stringsAsFactors = FALSE, header=T, sep = "\t") 
#File: picta_female_male_cov.txt

datacov_pic$MFlog2avg <- log2((datacov_pic$male)/(datacov_pic$female))
datasorted_pic <- datacov_pic[order(datacov_pic$chromosome,datacov_pic$start),]
autocov_pic <- datasorted_pic %>%
  filter_all(any_vars(datasorted_pic$chromosome != "LG12"))
autocov_pic <- datasorted_pic %>%
  filter(str_detect(datasorted_pic$chromosome, "^LG"))

autocov_pic <- na.omit(autocov_pic)

#calculate confidence interval for moving average for parae

MFautopermute_pic <- replicate(1000,mean(sample(autocov_pic$MFlog2avg,windowsizecov,replace = F)))
MFautoI25cov_pic <- quantile(MFautopermute_pic, c(.025, .5, .975))[[1]]
MFautoI25cov_pic
#[1] -0.1468106

MFautoI975cov_pic <- quantile(MFautopermute_pic, c(.025, .5, .975))[[3]]
MFautoI975cov_pic
#[1] 0.05038029


datacov_pic12 <- datasorted_pic %>%
  filter_all(any_vars(datasorted_pic$chromosome == "LG12"))


#LG12 also has been inverted and needs to be flipped:

pic_data12_edit <- datacov_pic12
pic_data12_edit$inv_start <- 32950001
pic_data12_edit$correct_start <- pic_data12_edit$inv_start - pic_data12_edit$start
pic_data12_edit$binstart <- pic_data12_edit$correct_start
datacov_pic12 <- pic_data12_edit


smoothlinefccov_pic = movingaverage(datacov_pic12$MFlog2avg, windowsizecov)

paneldatacov_pic <- as.data.frame(smoothlinefccov_pic)
paneldatacov_pic$startmb <- datacov_pic12$binstart
paneldatacov_pic <- na.omit(paneldatacov_pic)

ggplot(datacov_pic12, aes(x= binstart/1000000, y= MFlog2avg)) +
  annotate("rect", xmax = 34, xmin = 31, ymax = 1.2, ymin = -1.2, fill = "forestgreen", alpha = 0.3) +
  geom_rect(xmax=35,xmin=-10,ymax=MFautoI975cov_pic,ymin=MFautoI25cov_pic,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
  geom_line(data = paneldatacov_pic, aes(x=  paneldatacov_pic$startmb/1000000, y= paneldatacov_pic$smoothlinefccov), colour = "black", size=1.5) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,35)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +
  ylim(-1.2,1.2) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Coverage'))



######################################
############# P. parae ###############
######################################


datacov_par <- read.table(file.choose(), stringsAsFactors = FALSE, header=T, sep = "\t")
#File: female_male_parae_windows.txt

mean(datacov_par$fem_avg)
#[1] 14.36411
mean(datacov_par$male_avg)
#[1] 19.13979
### Difference between male_avg to female_avg is 1.332473087438066
### Adjust for sequence coverage differences between males and females
datacov_par$male_norm <-(datacov_par$male_avg/1.332473087438066)

datacov_par$MFlog2avg <- log2((datacov_par$male_norm)/(datacov_par$fem_avg))


datasorted_par <- datacov_par[order(datacov_par$CHROM,datacov_par$START),]
autocov_par <- datasorted_par %>%
  filter_all(any_vars(datasorted_par$CHROM != "LG12"))
autocov_par <- subset(autocov_par,!grepl("^Scaffold", autocov_par$CHROM))
autocov_par <- autocov_par %>%
  filter_all(any_vars(autocov_par$CHROM != "LG12"))


autocov_par <- na.omit(autocov_par)

#calculate confidence interval for moving average for parae

MFautopermute_par <- replicate(1000,mean(sample(autocov_par$MFlog2avg,windowsizecov,replace = T)))
MFautoI25cov_par <- quantile(MFautopermute_par, c(.025, .5, .975))[[1]]
MFautoI25cov_par
#[1] -0.07796683


MFautoI975cov_par <- quantile(MFautopermute_par, c(.025, .5, .975))[[3]]
MFautoI975cov_par
#[1] 0.03674249


datacov_par12 <- datasorted_par %>%
  filter_all(any_vars(datasorted_par$CHROM == "LG12"))


#LG12 also has been inverted and needs to be flipped:

par_data12_edit <- datacov_par12

par_data12_edit$inv_start <- 32950001
par_data12_edit$correct_start <- par_data12_edit$inv_start - par_data12_edit$START
par_data12_edit$binstart <- par_data12_edit$correct_start
datacov_par12 <- par_data12_edit
range(datacov_par12$MFlog2avg)
#[1] -0.973999  1.358927


smoothlinefccov_par = movingaverage(datacov_par12$MFlog2avg, windowsizecov)

paneldatacov_par <- as.data.frame(smoothlinefccov_par)
paneldatacov_par$startmb <- datacov_par12$binstart
paneldatacov_par <- na.omit(paneldatacov_par)

ggplot(datacov_par12, aes(x= binstart/1000000, y= MFlog2avg)) +
  annotate("rect", xmax = 35, xmin = 30, ymax = 1.2, ymin = -1.2, fill = "forestgreen", alpha = 0.3) +
  geom_rect(xmax=35,xmin=-10,ymax=MFautoI975cov_par,ymin=MFautoI25cov_par,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
  geom_line(data = paneldatacov_par, aes(x=  paneldatacov_par$startmb/1000000, y= paneldatacov_par$smoothlinefccov), colour = "black", size=1.5) +
  ylim(-1.2,1.2) +  coord_cartesian(xlim=c(0,35)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +  scale_x_continuous(breaks=seq(0,35,5)) +  xlab('Chromosome 12 (Mb)') +  ylab(expression('M:F log'[2]*' Coverage'))



