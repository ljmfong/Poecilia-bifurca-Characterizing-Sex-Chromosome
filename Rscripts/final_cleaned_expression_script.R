######################################
###### M:F Expression Analysis #######
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
movingaverage <- function(x, window) {
  roll_mean(x, n = window, fill = NA)  # Ensure correct arguments
}

windowsize = 40


#get data
data = read.table(file.choose(),header=T,stringsAsFactor=F,sep=",")
#File: gene_expression_position.txt

dim(data)
#[1] 16640     5
names(data)
#[1] "scaffold"   "gene"       "LG"         "pos"        "expression"


mean(data$expression)
#[1] 0.006657918
median(data$expression)
#[1] -0.03602223

#sort data
datasorted <- data[order(data$LG,data$pos),]
dim(datasorted)
#[1] 16640     5

#calculate confidence interval
dataminus12 <- datasorted[datasorted$LG !="LG12",]
dim(dataminus12)
#[1] 15962     5


MFpermutes <- replicate(1000,mean(sample(dataminus12$expression,windowsize,replace = TRUE)))
MFCI25 <- quantile(MFpermutes,c(.025, .5, .975))[[1]]
MFCI25
#[1] -0.5142813

MFCI975 <- quantile(MFpermutes,c(.025, .5, .975))[[3]]
MFCI975
#[1] 0.5988314

#Plot out for every other chromosome
data$LG <- factor(data$LG, paste0("LG",1:23), paste0("LG",1:23))


for (i in levels(data$LG)){
  chr_tester <- subset(data, LG == i)
  assign(paste0(i), chr_tester)
}

# List of objects starting with "LG"
data_list <- ls(pattern = "^LG")

#LG12 also has been inverted and needs to be flipped:

data12_edit <- LG12
data12_edit$inv_start <- 34233588
data12_edit$correct_start <- data12_edit$inv_start - data12_edit$pos
data12_edit$position_in_genome <- data12_edit$correct_start
data12_edit <- data12_edit[1:7]
LG12 <- data12_edit


generate_plot <- function(obj) {
  # Perform the operations
  smoothlinefc = movingaverage(obj$expression,windowsize)

  
  paneldata = as.data.frame(smoothlinefc)
  paneldata$startmb = obj$pos

  paneldata <- na.omit(paneldata)
  
  # Create the plot using ggplot
  plot <- ggplot(obj, aes(x = pos/1000000, y = expression)) +
    geom_rect(xmax=60,xmin=-10,ymax=MFCI975,ymin=MFCI25,fill="grey70", alpha=0.08)+
    geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1) +
    geom_line(data = paneldata, aes(x= startmb/1000000, y = smoothlinefc), size=1.5) +
    coord_cartesian(ylim=c(-0.5, 0.5)) +
    xlab(NULL) +
    ylab(NULL) +
    labs(title = obj$LG) +
    theme_classic() +
    theme(plot.title = element_text(size = 8, face = 'bold', hjust = 0.5),
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

final_plot <- grid.arrange(patchworkGrob(multi_panel_plot), left = "M:F log2 Expression (TPM)",
                           bottom = "Chromosome Length (Mb)")

final_plot


### Plotting the sex chromosome

smoothlinefccov = movingaverage(LG12$expression, windowsize)
smoothlinefemcov = movingaverage(LG12$expression, windowsize)
smoothlinemalcov = movingaverage(LG12$expression, windowsize)


paneldatacov <- as.data.frame(smoothlinefccov)
paneldatacov$startmb <- LG12$correct_start
names(paneldatacov)
# "smoothlinefccov" "startmb"
paneldatacov <- na.omit(paneldatacov)

# plots (you might have to change the axes to view your plot better - depends on the values plotted)
require(grid)

plotScaffold14 <- ggplot(LG12, aes(x= pos/1000000, y= expression)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFCI975,ymin=MFCI25,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1) +
  geom_line(data = paneldatacov, aes(x=  paneldatacov$startmb/1000000, y= paneldatacov$smoothlinefccov), colour = "black", size = 1.4) +
  coord_cartesian(xlim=c(0,35)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.8),
    axis.line.x = element_line(color="black", size = 0.8),
    axis.text.x = element_text(size=30),
    axis.text.y = element_text(size=30), axis.title = element_text(size = 30)
  ) +
  ylim(-2,2) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Expression (RPKM)'))

plotScaffold14


######################################
############# P. picta ###############
######################################

data_pic <- read.table(file.choose(),header=T,stringsAsFactor=F,sep=",")
#File: /home/ljmfong/bifurca_project/files_for_fig1/picta_alignment/rna/picta_all_individs_outfile_rpkm_fc_pos.txt


data_pic$expression <- log2((data_pic$avg_male_expr+0.01)/(data_pic$avg_female_expr+0.01))
data_pic$expression <- (data_pic$expression/1260333)*1e6
datasorted_pic <- data_pic[order(data_pic$chromosome,data_pic$position_in_genome),]
dataminus12_pic <- datasorted_pic[datasorted_pic$chromosome !="LG12",]
dim(dataminus12_pic)
#[1] 37839    14

MFpermutes_pic <- replicate(1000,mean(sample(dataminus12_pic$expression,windowsize,replace = TRUE)))
MFCI25_pic <- quantile(MFpermutes_pic,c(.025, .5, .975))[[1]]
MFCI25_pic
#[1] -0.1834453

MFCI975_pic <- quantile(MFpermutes_pic,c(.025, .5, .975))[[3]]
MFCI975_pic
#[1] 0.1351801

pic_LG12 <- subset(datasorted_pic, chromosome == "LG12")

#LG12 also has been inverted and needs to be flipped:

pic12_edit <- pic_LG12
pic12_edit$inv_start <- 32950001
pic12_edit$correct_start <- pic12_edit$inv_start - pic12_edit$position_in_genome
pic12_edit$position_in_genome <- pic12_edit$correct_start
pic_LG12 <- pic12_edit

smoothlinefccov_pic = movingaverage(pic_LG12$expression, windowsize)

paneldatacov_pic <- as.data.frame(smoothlinefccov_pic)
paneldatacov_pic$startmb <- pic_LG12$position_in_genome
paneldatacov_pic <- na.omit(paneldatacov_pic)

ggplot(pic_LG12, aes(x= position_in_genome/1000000, y= expression)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFCI975_pic,ymin=MFCI25_pic,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1) +
  geom_line(data = paneldatacov_pic, aes(x=  paneldatacov_pic$startmb/1000000, y= paneldatacov_pic$smoothlinefccov), colour = "black", size=1.2) +
  #coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,34)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +
  ylim(-0.5,0.5) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Expression (R{KM})'))


######################################
############# P. parae ###############
######################################

data_par <- read.table(file.choose(),header=T,stringsAsFactor=F,sep=",")
#File: parae_all_individs_outfile_rpkm_fc_pos.txt


datasorted_par <- data_par[order(data_par$chromosome,data_par$position_in_genome),]
dataminus12_par <- datasorted_par[datasorted_par$chromosome !="LG12",]
dim(dataminus12_par)
#[1] 26503     9


MFpermutes_par <- replicate(1000,mean(sample(dataminus12_par$expression,windowsize,replace = TRUE)))
MFCI25_par <- quantile(MFpermutes_par,c(.025, .5, .975))[[1]]
MFCI25_par
#[1] -0.0816574

MFCI975_par <- quantile(MFpermutes_par,c(.025, .5, .975))[[3]]
MFCI975_par
#[1] 0.05597327

par_LG12 <- subset(datasorted_par, chromosome == "LG12")

#LG12 also has been inverted and needs to be flipped:

par12_edit <- par_LG12
par12_edit$inv_start <- 34233588
par12_edit$correct_start <- par12_edit$inv_start - par12_edit$position_in_genome
par12_edit$position_in_genome <- par12_edit$correct_start
par_LG12 <- par12_edit

smoothlinefccov_par = movingaverage(par_LG12$expression, windowsize)

paneldatacov_par <- as.data.frame(smoothlinefccov_par)
paneldatacov_par$startmb <- par_LG12$position_in_genome
paneldatacov_par <- na.omit(paneldatacov_par)

ggplot(par_LG12, aes(x= position_in_genome/1e6, y= expression)) +
  geom_rect(xmax=60,xmin=0,ymax=MFCI975_par,ymin=MFCI25_par,fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1) +
  geom_line(data = paneldatacov_par, aes(x=  paneldatacov_par$startmb/1e6, y= paneldatacov_par$smoothlinefccov), colour = "black", size=1.2) +
  coord_cartesian(xlim=c(0,35)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.8),
    axis.line.x = element_line(color="black", size = 0.8),
    axis.text.x = element_text(size=30),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 30)
  ) +
  ylim(-0.5,0.5) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Expression (RPKM)'))



