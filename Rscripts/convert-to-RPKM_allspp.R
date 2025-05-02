
#Convert counts to RPKM

rm(list=ls())
ls() 


library(edgeR)
library(tibble)

data <- read.table(file.choose(),stringsAsFactors=F,header=T, row.names=1)
#data <- read.table("all_bifurca_read_count.txt",stringsAsFactors=F,header=T, row.names=1)
#data FORMAT: Geneid	Sample1	Sample2
#data FORMAT: MSTRG.43175	145.0	285.0
dim(data)
#[1] 27166     6
names(data)

conditions <- factor(c("F", "F", "F", "M", "M", "M"))

#Extract RPKM
expr <- DGEList(counts=data)
gene_length <- read.table(file.choose(),stringsAsFactors=F)
#gene_length <- read.table("bifurca_gene_length.txt",stringsAsFactors=F)

names(gene_length)
#[1] "V1" "V2"

dim(gene_length)
#[1] 27166     2


expressed_genes <- rownames(data)

gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE

bifurca_all_rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector) #This will normalize your genes for gene length
str(bifurca_all_rpkm)

write.table(bifurca_all_rpkm, file="bifurca_all_rpkm.txt",quote=F, sep="\t")

#Filter genes for low expression: 2rpkm
#Done using the python script: 15.filter-expression-2rpkm-halformoreonesex.py


##################
#### P. parae ####
##################

data <- read.table(file.choose(),stringsAsFactors=F,header=T, row.names=1)
#data <- read.table("parae_read_counts.txt",stringsAsFactors=F,header=T, row.names=1)
#data FORMAT: Geneid	Sample1	Sample2
#data FORMAT: MSTRG.43175	145.0	285.0
dim(data)
#[1] 401776     6
names(data)

conditions <- factor(c("F", "F", "F", "M", "M", "M"))

#Extract RPKM
expr <- DGEList(counts=data)
gene_length <- read.table(file.choose(),stringsAsFactors=F)
#Use the picta BED file and convert to gene lengths
#gene_length <- read.table("Poecilia_picta_gene_length.txt",stringsAsFactors=F)

names(gene_length)
#[1] "V1" "V2"

dim(gene_length)
#[1] 27764     2


expressed_genes <- rownames(data)

gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE

parae_all_rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector) #This will normalize your genes for gene length
str(parae_all_rpkm)

getwd()

write.table(parae_all_rpkm, file="parae_all_rpkm.txt",quote=F, sep="\t")
#Move the file to its respective folder (moved to count_extraction_parae)

#Do the same thing but for picta gene counts:
data <- read.table(file.choose(),stringsAsFactors=F,header=T, row.names=1)
#picta_extract_all_read_count.txt
conditions <- factor(c("F", "F", "F", "M", "M", "M"))
expr <- DGEList(counts=data)
expressed_genes <- rownames(data)

gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))


picta_all_rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector) #This will normalize your genes for gene length
str(picta_all_rpkm)

getwd()

write.table(picta_all_rpkm, file="picta_all_rpkm.txt",quote=F, sep="\t")
#Move to the respective folder

#Filter genes for low expression: 2rpkm
#Done using the python script: 15.filter-expression-2rpkm-halformoreonesex.py
