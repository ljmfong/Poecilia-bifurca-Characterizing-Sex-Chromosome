###Unsplit fdrs

data = read.table(file.choose(),stringsAsFactors=F,header=F,sep="\t")
#"bifurca_project/continue_ase/female_pvalues.txt

head(data)

p = data$V2
head(p)


p_fdr = p.adjust(p, method=c("fdr"), n= length(p))
head(p_fdr)
head(p)

data$pvaladj = p_fdr
head(data)

getwd()


write.table(data, file="females_list_pvaluesadj.txt",sep="\t")

#Move this file into the respective folder


data = read.table(file.choose(),stringsAsFactors=F,header=F,sep="\t")
#"bifurca_project/continue_ase/male_pvalues.txt

head(data)

p = data$V2
head(p)


p_fdr = p.adjust(p, method=c("fdr"), n= length(p))
head(p_fdr)
head(p)

data$pvaladj = p_fdr
head(data)

getwd()


write.table(data, file="males_list_pvaluesadj.txt",sep="\t")
#Move this file into the respective folder
