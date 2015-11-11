setwd("Z:/Chirayu/My Projects/DREAM7/Validation Set/SECONDARY (STATISTICAL) ANALYSIS/")
ds<-read.csv(file="Sig genes low ic20.txt", sep="\t")

drugs<-c("DEME","INDO","MEBE","MICO","NISO","SULC","TIOC","TORE","DIGI","DYCL","MITO","TENI","TOPO")


#SYNERGISM
zz<-file("results_syn.txt")
for ( i in 1:length(drugs)){
	drug1<-as.character(drugs[i])
	drug1pvalcol<-as.character(paste(drug1, "_PVAL", sep=""))
	no_sig_drg1<-length(which(ds[,drug1pvalcol]<=0.05))
	drug1fccol<-as.character(paste(drug1,"_FC", sep=""))
	
	for (j in 1:length(drugs)){
		if (j > i){
			drug2<-as.character(drugs[j])
			drug2pvalcol<-as.character(paste(drug2,"_PVAL", sep=""))
			no_sig_drg2<-length(which(ds[,drug2pvalcol]<=0.05))
			drug2fccol<-as.character(paste(drug2,"_FC", sep=""))
			total=0
			for ( m in 1:nrow(ds)){
				 if(sign(ds[drug1fccol][m,])==sign(ds[drug2fccol][m,]) && abs(ds[drug1fccol][m,])>1.1 && abs(ds[drug2fccol][m,])>1.1)
				 {
				 total=total+1
				 #print(total)
				 }
			 }
		print(paste(drug1," - ",drug2,": ",total," - ", (no_sig_drg1+no_sig_drg2), sep=""))
		writeLines(paste(drug1,drug2,total, sep="\t"), zz)
		}
	}
}

close(zz)



#ANTAGONISM 
zz<-file("results_anta.txt")
for ( i in 1:length(drugs)){
	drug1<-as.character(drugs[i])
	drug1pvalcol<-as.character(paste(drug1, "_PVAL", sep=""))
	no_sig_drg1<-length(which(ds[,drug1pvalcol]<=0.05))
	drug1fccol<-as.character(paste(drug1,"_FC", sep=""))
	
	for (j in 1:length(drugs)){
		if (j > i){
			drug2<-as.character(drugs[j])
			drug2pvalcol<-as.character(paste(drug2,"_PVAL", sep=""))
			no_sig_drg2<-length(which(ds[,drug2pvalcol]<=0.05))
			drug2fccol<-as.character(paste(drug2,"_FC", sep=""))
			denom<-no_sig_drg1+no_sig_drg2
			total=0
			for ( m in 1:nrow(ds)){
				 if(sign(ds[drug1fccol][m,])!=sign(ds[drug2fccol][m,]) && abs(ds[drug1fccol][m,])>1.1 && abs(ds[drug2fccol][m,])>1.1)
				 {
				 total=total+1
				 #print(total)
				 }
			 }
		
		#print(paste(drug1," - ",drug2,": ",total, sep=""))
		print(paste(drug1," - ",drug2,": ",total," - ", (no_sig_drg1+no_sig_drg2), sep=""))
		#writeLines(paste(drug1,drug2,total, sep="\t"), zz)
		}
	}
}
close(zz)



#ADDITIVISM
zz<-file("results_add.txt")
for ( i in 1:length(drugs)){
	drug1<-as.character(drugs[i])
	drug1pvalcol<-as.character(paste(drug1, "_PVAL", sep=""))
	drug1fccol<-as.character(paste(drug1,"_FC", sep=""))
	
	for (j in 1:length(drugs)){
		if (j > i){
			drug2<-as.character(drugs[j])
			drug2pvalcol<-as.character(paste(drug2,"_PVAL", sep=""))
			drug2fccol<-as.character(paste(drug2,"_FC", sep=""))
			total=0
			for ( m in 1:nrow(ds)){
				 if(sign(ds[drug1fccol][m,])==sign(ds[drug2fccol][m,]) && (ds[drug1pvalcol][m,]<=0.05 | ds[drug2pvalcol][m,]<=0.05))
				 {
				 total=total+1
				 #print(total)
				 }
			 }
		
		print(paste(drug1," - ",drug2,": ",total, sep=""))
		#writeLines(paste(drug1,drug2,total, sep="\t"), zz)
		}
		
	}
	
}

close(zz)
