#cov_fc_mov_av

rm(list=ls())
ls() 

library(ggplot2)
library(reshape)
library(gridExtra)
library(dplyr)
library(ggthemes)
library(extrafont)
loadfonts()
library(RcppRoll)
library(zoo)
library(ggsave)

#define a moving average function
movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill = NA)
}

windowsizecov = 50

#get data
datacov <- read.csv(file.choose(), stringsAsFactors = FALSE, header=T) 
#This should also have positional info 
#In this case, grab the start and end of each Scaffold using Dave's picta genome
#It is worth noting that in the RagTag assembly, there are more scaffolds than what were put together in Dave's P. picta genome
#Therefore, to try and make a comparison with only the scaffolds that were assembled
#Get only the values from Scaffold 1 - Scaffold 23 --> Manually done in the excel file
#File is \\files.zoology.ubc.ca\ljmfong\flex\bifurca_project\coverage_analysis\MF_log_avg_cov_pos_23scaffs.csv

colnames(datacov)[1] <- "Scaffold"

datacov$Mlogaverage <- log(datacov$Mcovavg)
datacov$Flogaverage <- log(datacov$Fcovavg)
#datacov$MFLogaverage <- datacov$Mlogaverage - datacov$Flogaverage
datacov$MFlog2avg <- log2((datacov$Mcovavg)/(datacov$Fcovavg))

autocov <- datacov %>%
  filter_all(any_vars(datacov$Scaffold != "Scaffold_14_RagTag"))

#We are assuming that Scaffold 14 (which is LG12 in P. picta is the sex chromosome)

  
dim(datacov)
#[1] 16518  8
dim(autocov)
#[1] 15833  8

names(datacov)
# [1] "Scaffold"  "binstart"  "binend"    "Fcovavg"   "Mcovavg"   "MFlog2avg"
#[7] "Mlogaverage" "Flogaverage"

mean(autocov$MFlog2avg)
# [1] 0.0917997


#sort data - in this case, the chromosome is already sorted into the bins of 50kB
autocovsorted <- autocov[order(autocov$Scaffold,autocov$binstart),]
dim(autocovsorted)
# [1] 15833     8

#calculate confidence interval for moving average

MFautopermute <- replicate(1000,mean(sample(autocovsorted$MFlog2avg,windowsizecov,replace = FALSE)))
MFautoI25cov <- quantile(MFautopermute, c(.025, .5, .975))[[1]]
MFautoI25cov
#[1] 0.0236532
MFautoI975cov <- quantile(MFautopermute, c(.025, .5, .975))[[3]]
MFautoI975cov
#[1] 0.1421708

#Plotting Scaffold_14/LG12

datacovScaffold14 <- datacov %>%
  filter_all(any_vars(datacov$Scaffold == "Scaffold_14_RagTag"))
  
datacovScaffold14 <- datacovScaffold14[order(datacovScaffold14$binstart,as.numeric(levels(datacovScaffold14$binstart))[datacovScaffold14$binstart]),]

dim(datacovScaffold14)
# [1] 685   8

mean(datacovScaffold14$MFlog2avg)
# [1] -0.7815422


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
  geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacov, aes(x=  paneldatacov$startmb/1000000, y= paneldatacov$smoothlinefccov), colour = "black", size=0.6) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,32)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  ylim(-2,2) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('LG12 Start position (Mb)') +
  ylab(expression('M:F log'[2]*' coverage'))


plotScaffold14

plotsexcov <- ggplot(datacovScaffold14, aes(x= binstart/1000000, y= Flogaverage)) + 
  geom_point(colour="#e31a1c",fill="#e31a1c", alpha=0.5, cex=0.4) +
  geom_point(data= datacovScaffold14, aes(x= binstart/1000000, y= Mlogaverage),colour="#1f78b4",fill="#1f78b4", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$startmb/1000000, y=paneldatacov$females), colour="#e31a1c",, size=0.6) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$startmb/1000000, y=paneldatacov$males), colour="#1f78b4",, size=0.6) +
  coord_cartesian(ylim=c(0,4)) +
  coord_cartesian(xlim=c(0,32)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('LG12 Start position (Mb)') +
  ylab(expression('Log'[2]*' coverage'))

plotsexcov

require(cowplot)
plot_grid(plotScaffold14, plotsexcov, labels = c('A', 'B'),ncol = 1,nrow=2, align="v")


#Orient the sex chromosome into the correct orientation

paneldatacov$inv_value <- 32950001
paneldatacov$inv_start <- paneldatacov$inv_value - paneldatacov$startmb
datacovScaffold14$inv_value <- 32950001
datacovScaffold14$inv_start <- datacovScaffold14$inv_value - datacovScaffold14$binstart



plotScaffold14 <- ggplot(datacovScaffold14, aes(x= inv_start/1000000, y= MFlog2avg)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacov, aes(x=  paneldatacov$inv_start/1000000, y= paneldatacov$smoothlinefccov), colour = "black", size=0.6) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,32)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  ylim(-2,2) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('LG12 Start position (Mb)') +
  ylab(expression('M:F log'[2]*' coverage'))

plotScaffold14


plotsexcov <- ggplot(datacovScaffold14, aes(x= inv_start/1000000, y= Flogaverage)) + 
  geom_point(colour="#e31a1c",fill="#e31a1c", alpha=0.5, cex=0.4) +
  geom_point(data= datacovScaffold14, aes(x= inv_start/1000000, y= Mlogaverage),colour="#1f78b4",fill="#1f78b4", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$inv_start/1000000, y=paneldatacov$females), colour="#e31a1c",, size=0.6) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$inv_start/1000000, y=paneldatacov$males), colour="#1f78b4",, size=0.6) +
  coord_cartesian(ylim=c(0,4)) +
  coord_cartesian(xlim=c(0,32)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('LG12 Start position (Mb)') +
  ylab(expression('Log'[2]*' coverage'))

plotsexcov

plot_grid(plotScaffold14, plotsexcov, labels = c('A', 'B'),ncol = 1,nrow=2, align="v")


#Test an autosome to compare:

Scaff1 <- datacov %>%
  filter_all(any_vars(datacov$Scaffold == "Scaffold_1_RagTag"))

Scaff1 <- Scaff1[order(Scaff1$binstart,as.numeric(levels(Scaff1$binstart))[Scaff1$binstart]),]

mean(Scaff1$MFlog2avg)
# [1] 0.09510463

smoothlinefccovScaff1 = movingaverage(Scaff1$MFlog2avg, windowsizecov)
smoothlinefemcovScaff1 = movingaverage(Scaff1$Flogaverage, windowsizecov)
smoothlinemalcovScaff1 = movingaverage(Scaff1$Mlogaverage, windowsizecov)

paneldatacovScaff1 <- as.data.frame(smoothlinefccovScaff1)
paneldatacovScaff1$startmb <- Scaff1$binstart
paneldatacovScaff1$females <- smoothlinefemcovScaff1
paneldatacovScaff1$males <- smoothlinemalcovScaff1
paneldatacovScaff1 <- na.omit(paneldatacovScaff1)

plotScaffold1 <- ggplot(Scaff1, aes(x= binstart/1000000, y= MFlog2avg)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacovScaff1, aes(x= paneldatacovScaff1$startmb/1000000, y= paneldatacovScaff1$smoothlinefccovScaff1), colour = "black", size=0.6) +
  coord_cartesian(ylim=c(-1,1)) +
  coord_cartesian(xlim=c(0,50)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  ylim(-1,1) + 
  scale_x_continuous(breaks=seq(0,50,5)) +
  xlab('LG2 Start position (Mb)') +
  ylab(expression('M:F log'[2]*' coverage'))

plotScaffold1

plotsexcovScaff1 <- ggplot(Scaff1, aes(x= binstart/1000000, y= Flogaverage)) + 
  geom_point(colour="#e31a1c",fill="#e31a1c", alpha=0.5, cex=0.4) +
  geom_point(data= Scaff1, aes(x= binstart/1000000, y= Mlogaverage),colour="#1f78b4",fill="#1f78b4", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacovScaff1, aes(x= paneldatacovScaff1$startmb/1000000, y=paneldatacovScaff1$females), colour="#e31a1c",, size=0.6) +
  geom_line(data = paneldatacovScaff1, aes(x= paneldatacovScaff1$startmb/1000000, y=paneldatacovScaff1$males), colour="#1f78b4",, size=0.6) +
  coord_cartesian(ylim=c(0,4)) +
  coord_cartesian(xlim=c(0,50)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,50,5)) +
  xlab('LG2 Start position (Mb)') +
  ylab(expression('Log'[2]*' coverage'))

plotsexcovScaff1

plot_grid(plotScaffold1, plotsexcovScaff1, labels = c('A', 'B'),ncol = 1,nrow=2, align="v")


#Test a smaller autosome - LG23/SCaffold 23

Scaff23 <- datacov %>%
  filter_all(any_vars(datacov$Scaffold == "Scaffold_23_RagTag"))

Scaff23 <- Scaff23[order(Scaff23$binstart,as.numeric(levels(Scaff23$binstart))[Scaff23$binstart]),]

mean(Scaff23$MFlog2avg)
# [1] 0.09118017

smoothlinefccovScaff23 = movingaverage(Scaff23$MFlog2avg, windowsizecov)
smoothlinefemcovScaff23 = movingaverage(Scaff23$Flogaverage, windowsizecov)
smoothlinemalcovScaff23 = movingaverage(Scaff23$Mlogaverage, windowsizecov)

paneldatacovScaff23 <- as.data.frame(smoothlinefccovScaff23)
paneldatacovScaff23$startmb <- Scaff23$binstart
paneldatacovScaff23$females <- smoothlinefemcovScaff23
paneldatacovScaff23$males <- smoothlinemalcovScaff23
paneldatacovScaff23 <- na.omit(paneldatacovScaff23)

plotScaffold23 <- ggplot(Scaff23, aes(x= binstart/1000000, y= MFlog2avg)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacovScaff23, aes(x= paneldatacovScaff23$startmb/1000000, y= paneldatacovScaff23$smoothlinefccovScaff23), colour = "black", size=0.6) +
  coord_cartesian(ylim=c(-1,1)) +
  coord_cartesian(xlim=c(0,25)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  ylim(-1,1) + 
  scale_x_continuous(breaks=seq(0,25,5)) +
  xlab('LG23 Start position (Mb)') +
  ylab(expression('M:F log'[2]*' coverage'))

plotScaffold23

plotsexcovScaff23 <- ggplot(Scaff23, aes(x= binstart/1000000, y= Flogaverage)) + 
  geom_point(colour="#e31a1c",fill="#e31a1c", alpha=0.5, cex=0.4) +
  geom_point(data= Scaff23, aes(x= binstart/1000000, y= Mlogaverage),colour="#1f78b4",fill="#1f78b4", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacovScaff23, aes(x= paneldatacovScaff23$startmb/1000000, y=paneldatacovScaff23$females), colour="#e31a1c",, size=0.6) +
  geom_line(data = paneldatacovScaff23, aes(x= paneldatacovScaff23$startmb/1000000, y=paneldatacovScaff23$males), colour="#1f78b4",, size=0.6) +
  coord_cartesian(ylim=c(0,4)) +
  coord_cartesian(xlim=c(0,25)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,25,5)) +
  xlab('LG23 Start position (Mb)') +
  ylab(expression('Log'[2]*' coverage'))

plotsexcovScaff23

plot_grid(plotScaffold23, plotsexcovScaff23, labels = c('A', 'B'),ncol = 1,nrow=2, align="v")




