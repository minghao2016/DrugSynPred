setwd("Z:/Chirayu/My Projects/DREAM7/Validation Set/SECONDARY (STATISTICAL) ANALYSIS/")
ds<-read.csv(file="Differentially Expressed Genes.txt", sep="\t")

low<-c("DIGI","DYCL","MITO","TENI","TOPO")
high<-c("DEME","INDO","MEBE","MICO","NISO","SULC","TIOC","TORE")

# genes_high<-c()
# for (j in 1:length(high)){
	# name<-paste(as.character(high[j]), "_PVAL", sep="")
	# genes_high<-append(genes_high, as.character(ds[,1][which(ds[,name] <=0.05)]))
# }
# genes_high<-unique(genes_high)

genes_low<-c()
for (j in 1:length(low)){
	name<-paste(as.character(low[j]), "_PVAL", sep="")
	genes_low<-append(genes_low, as.character(ds[,1][which(ds[,name] <=0.05)]))
}
genes_low<-unique(genes_low)

length(genes_low)
write.table(as.data.frame(genes_low), file="SIG genes p val 5perc.txt", sep="\t", row.names=F, quote=F)
