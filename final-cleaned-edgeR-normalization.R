library("edgeR")

data <- read.table(file.choose(),stringsAsFactors=F,header=T, row.names=1)
#File: extract_all_read_count.txt

#data FORMAT: Geneid	Sample1	Sample2
#data FORMAT: MSTRG.43175	145.0	285.0
names(data)
#[1] "bif_female_1_sorted" "bif_female_2_sorted" "bif_female_3_sorted" "bif_male_1_sorted"  
#[5] "bif_male_2_sorted"   "bif_male_3_sorted"  
dim(data)
#[1] 27166     6
conditions <- factor(c("F","F","F","M","M","M"))
#conditions eg conditions <- factor(c("M","M","M","M","M","M","M","M","M","M","F","F","F","F","M"))
length(conditions)
#[1] 6

#Check raw read count data
expr <- DGEList(counts=data,group=conditions)
group = gl(6,1,labels=c("Female1","Female2","Female3","Male1","Male2","Male3"))
plotMDS(expr,labels=group)
plotSmear(expr,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
abline(h=0,col="gray21",lty=2)
abline(h=2,col="gray",lty=2)
abline(h=-2,col="gray",lty=2)

cpm_expr <- cpm(expr)
sample1 <- density(log2(cpm_expr[,1]))
sample2 <- density(log2(cpm_expr[,2]))
sample3 <- density(log2(cpm_expr[,3]))
sample4 <- density(log2(cpm_expr[,4]))
sample5 <- density(log2(cpm_expr[,5]))
sample6 <- density(log2(cpm_expr[,6]))
plot(sample1, xlab = "CPM", ylab = "Density",legend(7,0.20,c("Male_bifurca","Female_bifurca"),lty=c(1,1),lwd=c(2,2),col=c("blue","red")),ylim=c(0,0.22),type="l",lwd=2,main="Raw log2 cpm", col="red")
lines(sample2, type="l",lwd=2,col="red2")
lines(sample3, type="l",lwd=2,col="red3")
lines(sample4, type="l",lwd=2,col="blue1")
lines(sample5, type="l",lwd=2,col="blue2")
lines(sample6, type="l",lwd=2,col="blue3")

#Check normalised read count data
expr <- DGEList(counts=data,group=conditions)
norm_expr <- calcNormFactors(expr)
group = gl(6,1,labels=c("Female1","Female2","Female3","Male1","Male2","Male3"))
plotMDS(norm_expr,labels=group)
plotSmear(norm_expr,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
abline(h=0,col="gray21",lty=2)
abline(h=2,col="gray",lty=2)
abline(h=-2,col="gray",lty=2)

cpm_norm_expr <- cpm(norm_expr)
sample1 <- density(log2(cpm_norm_expr[,1]))
sample2 <- density(log2(cpm_norm_expr[,2]))
sample3 <- density(log2(cpm_norm_expr[,3]))
sample4 <- density(log2(cpm_norm_expr[,4]))
sample5 <- density(log2(cpm_norm_expr[,5]))
sample6 <- density(log2(cpm_norm_expr[,6]))

plot(sample1, xlab = "CPM", ylab = "Density",legend(7,0.20,c("Male_bifurca","Female_bifurca"),lty=c(1,1),lwd=c(2,2),col=c("blue","red")),ylim=c(0,0.22),type="l",lwd=2,main="Raw log2 cpm",col="red")
lines(sample2, type="l",lwd=2,col="red2")
lines(sample3, type="l",lwd=2,col="red3")
lines(sample4, type="l",lwd=2,col="blue")
lines(sample5, type="l",lwd=2,col="blue2")
lines(sample6, type="l",lwd=2,col="blue3")

#Normalise and extract rpkm
expr <- DGEList(counts=data)
norm_expr <- calcNormFactors(expr)
gene_length <- read.table(file.choose(),stringsAsFactors=F)
#File: /home/ljmfong/bifurca_project/07.gene_expression/count_extractions/bifurca_gene_length.txt
names(gene_length)
#[1] "V1" "V2"
dim(gene_length)
#[1] 27166     2
expressed_genes <- rownames(data)
length(expressed_genes)
#[1] 27166
gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#[1] TRUE
#should print TRUE
rpkm_norm <- rpkm(norm_expr, log=FALSE,gene.length=gene_length_vector)

getwd()

write.table(rpkm_norm, file="bifurca_read_counts_rpkm_normalised.txt",quote=F, sep="\t")
#log2 rpkm
log_rpkm_norm = log2(rpkm_norm+1)
write.table(log_rpkm_norm, file="bifurca_read_counts_rpkm_normalised_log2.txt",quote=F, sep=",")

#Move to appropriate directory: normalized_reads

#############################################
### Do the same for P. picta and P. parae ###
#############################################




