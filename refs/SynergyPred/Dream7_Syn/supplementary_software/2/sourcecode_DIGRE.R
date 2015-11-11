Title: The challenge of predicting synergistic and antagonistic compound-pair activity from individual compound perturbations

Mukesh Bansal, Jichen Yang, Charles Karan, Michael Patrick Menden, James C. Costello, Hao Tang, Guanghua Xiao, Yajuan Li, Jeffrey Allen, Rui Zhong, Beibei Che8, Minsoo Kim, Tao Wang, Laura Heiser, Ronald Realubit, Michela Mattioli, Mariano Alvarez, Yao Shen, NCI-DREAM community, Daniel Gallahan, Dinah Singer, Julio Saez-Rodriguez, Yang Xie, Gustavo Stolovitzky, Andrea Califano


####################################################################
# source/pseudo code for DIGRE Model        
# Questions or comments for the code, software and website, 
# please contact: jichenyang@gmail.com and bo.yao@utsouthestern.edu 
####################################################################


The main hypothesis for DIGRE is that the transcriptomic residual effects induced by treatment with a previous drug(s) modulate the effect of a successively administered drug(s), e.g., suppose a cell line is treated with compounds A and B sequentially. The compound A, beside its own drug effect, also alters the cell line's genomic context, which further modulates the cells's response to compound B. Therefore, it is important to estimate this residual effect, in order to predict the net effect of the drug combination.

############### The scoring function  ################
Function "scoring" includes three parts:
#Part 1, derive similarity score between drug A and B. 
#Part 2, estimate drug induced genomic residual effect by integrating similarity score and 
#drug response curves. Part 3, Predict the synergistic score. 
#Function "scoring" includes three inputs 'ICDiff', 'resPath1', and 'resPath2'. 'ICDiff' represents 
#the differences of gene expression (in log ratio) before
#and after drug treatment. 
#For example:
#         AffyID  Genename Aclacinomycin A Blebbistatin Camptothecin Cycloheximide
#   11715100_at  HIST1H3G     -0.13500000  -0.46166667  -0.12166667    0.09833333
# 11715103_x_at TNFAIP8L1      0.02750000  -0.11250000   0.08416667   -0.12916667
# 11715104_s_at     OTOP2     -0.04958333  -0.18625000  -0.01958333   -0.04958333
#   11715105_at  C17orf78     -0.01708333   0.03291667  -0.05375000   -0.09041667
#  Doxorubicin-hydrochloride   Etoposide Geldanamycin H-7, Dihydrochloride Methotrexate
#               -0.08833333 -0.18833333  -0.04833333          -0.12500000   0.16166667
#               -0.06583333  0.15416667   0.07416667           0.13750000   0.13083333
#                0.03041667  0.08375000  -0.16958333           0.06375000   0.04708333
#               -0.01041667 -0.02041667   0.11958333          -0.03041667   0.14958333
#  Mitomycin C   Monastrol   Rapamycin Trichostatin A Vincristine
# -0.21500000  0.13166667 -0.07833333    -0.01500000 -0.30166667
# -0.04250000 -0.03916667 -0.04250000     0.73750000 -0.06916667
# -0.27291667 -0.16291667  0.11041667    -0.03958333 -0.02291667
# -0.16041667 -0.14375000 -0.05041667    -0.17041667 -0.04041667
#
#'resPath1' and 'resPath2' are the folders to save intermediate results.
#In addition, the input includes 'bigGraph', 'bigGraph2'
#, 'bigGraphnodes', and 'bigGraphnodes2', which are pre-identified cell growth related 
#KEGG pathways (CGP), global pathways (GP), CGP's nodes, 
#and GP's nodes, respectively. The KEGG pathways 
#included in CGP and GP were explicitly explained
#in the method description. They can be easily derived by the R package: KEGGgraph. 


scoring<-function(ICDiff,resPath1,resPath2){
  #### preparation: remove redundant data, GeneID to KEGGID
  resIC20<-list()
  names3<-colnames(ICDiff)[3:length(ICDiff[1,])]
  for(i in 1:(length(ICDiff[1,])-2)){
    print(i)
    for(j in 1:(length(ICDiff[1,])-2)){
      if(i==j) next
      drugU1<-ICDiff[ICDiff[,i+2]>fold,c(2,i+2)]
      drugU2<-ICDiff[ICDiff[,j+2]>fold,c(2,j+2)]
      drugU1<-drugU1[!duplicated(drugU1[,1]),]
      drugU2<-drugU2[!duplicated(drugU2[,1]),]
      U1<-drugU1[,1];U2<-drugU2[,1]
      U1ID<-sapply(mget(U1,org.Hs.egSYMBOL2EG,ifnotfound=NA),"[[",1)
      U1ID<-U1ID[!is.na(U1ID)]
      U2ID<-sapply(mget(U2,org.Hs.egSYMBOL2EG,ifnotfound=NA),"[[",1)
      U2ID<-U2ID[!is.na(U2ID)]
      U1ID<-translateGeneID2KEGGID(U1ID, organism="hsa")
      U2ID<-translateGeneID2KEGGID(U2ID, organism="hsa")
      ## differentially expressed genes after treatment by Drug D1( i.e. drug A) or drug D2 (i.e. drug B)
      drugD1<-ICDiff[ICDiff[,i+2]<(-fold),c(2,i+2)]
      drugD2<-ICDiff[ICDiff[,j+2]<(-fold),c(2,j+2)]
      drugD1<-drugD1[!duplicated(drugD1[,1]),]
      drugD2<-drugD2[!duplicated(drugD2[,1]),]
      D1<-drugD1[,1];D2<-drugD2[,1]
      D1ID<-sapply(mget(D1,org.Hs.egSYMBOL2EG,ifnotfound=NA),"[[",1)
      D1ID<-D1ID[!is.na(D1ID)]
      D2ID<-sapply(mget(D2,org.Hs.egSYMBOL2EG,ifnotfound=NA),"[[",1)
      D2ID<-D2ID[!is.na(D2ID)]
      D1ID<-translateGeneID2KEGGID(D1ID, organism="hsa")
      D2ID<-translateGeneID2KEGGID(D2ID, organism="hsa")
      
      Ushared<-U1ID[U1ID %in% U2ID]
      U1ID<-U1ID[!U1ID %in% Ushared]
      U2ID<-U2ID[!U2ID %in% Ushared]
      Dshared<-D1ID[D1ID %in% D2ID]
      D1ID<-D1ID[!D1ID %in% Dshared]
      D2ID<-D2ID[!D2ID %in% Dshared]
      U1nodes<-U1ID[U1ID %in% bigGraphnodes]
      U2nodes<-U2ID[U2ID %in% bigGraphnodes2]
      D1nodes<-D1ID[D1ID %in% bigGraphnodes]
      D2nodes<-D2ID[D2ID %in% bigGraphnodes2]

      #### Part 1, derive similarity score between drug A and B
      
      ############### for up-regulated genes in drug A (URG_A) ####################
      U1graph<-subKEGGgraph(U1nodes,bigGraph)
      U1edges<-getKEGGedgeData(U1graph)
      U1p<-0;U1n<-0;
      if(length(U1edges)>0){
        for(k in 1:length(U1edges)){
          if(U1edges[[k]]@subtype[[1]]@name=="activation"
             | U1edges[[k]]@subtype[[1]]@name=="expression") U1p=U1p+1
          if(U1edges[[k]]@subtype[[1]]@name=="inhibition") U1n=U1n+1
        }
      }
      ####UP+UP
      UID<-c(U1ID, U2ID)
      Up<-0;Un<-0;
      sharedPositive<-Ushared[Ushared %in% bigGraphnodes]
      Unodes<-UID[UID %in% bigGraphnodes2]
      Ugraph<-subKEGGgraph(Unodes,bigGraph2)
      Uedges<-getKEGGedgeData(Ugraph)
      if(length(Uedges)>0){
        for(k in 1:length(Uedges)){
          if(strsplit(names(Uedges[k]),"~")[[1]][2] %in% U1nodes){
            if(Uedges[[k]]@subtype[[1]]@name=="activation"
               | Uedges[[k]]@subtype[[1]]@name=="expression") Up=Up+1
            if(Uedges[[k]]@subtype[[1]]@name=="inhibition") Un=Un+1
          }
        }
      }
      
      ####UP+DOWN
      UDID<-c(U1ID,D2ID)
      UDp<-0;UDn<-0;
      UDnodes<-UDID[UDID %in% bigGraphnodes]
      UDn<-sum(as.numeric(duplicated(UDnodes)))
      UDnodes<-UDID[UDID %in% bigGraphnodes2]
      common<-U1ID[U1ID %in% D2ID]
      UDnodes<-UDnodes[!UDnodes %in% common]
      UDnodes<-UDnodes[!duplicated(UDnodes)]
      UDgraph<-subKEGGgraph(UDnodes,bigGraph2)
      UDedges<-getKEGGedgeData(UDgraph)
      if(length(UDedges)>0){
        for(k in 1:length(UDedges)){
          if(strsplit(names(UDedges[k]),"~")[[1]][2] %in% U1nodes
             & strsplit(names(UDedges[k]),"~")[[1]][1] %in% D2ID){
            if(UDedges[[k]]@subtype[[1]]@name=="activation"
               | UDedges[[k]]@subtype[[1]]@name=="expression") UDn=UDn+1
            if(UDedges[[k]]@subtype[[1]]@name=="inhibition") UDp=UDp+1
          }
        }
      }
      
      Upositive=Up+UDp-U1p
      Unegtive=Un+UDn-U1n
      
      #### for DRG_A
      ############### for drug1 down-regulated ####################
      D1graph<-subKEGGgraph(D1nodes,bigGraph)
      D1edges<-getKEGGedgeData(D1graph)
      D1p<-0;D1n<-0;
      if(length(D1edges)>0){
        for(k in 1:length(D1edges)){
          if(D1edges[[k]]@subtype[[1]]@name=="activation"
             | D1edges[[k]]@subtype[[1]]@name=="expression") D1p=D1p+1
          if(D1edges[[k]]@subtype[[1]]@name=="inhibition") D1n=D1n+1
        }
      }
      ####DOWN+UP
      DUID<-c(D1ID, U2ID)
      DUp<-0;DUn<-0;
      DUnodes<-DUID[DUID %in% bigGraphnodes]
      DUn=DUn+sum(as.numeric(duplicated(DUnodes)))
      DUnodes<-DUID[DUID %in% bigGraphnodes2]
      common<-D1ID[D1ID %in% U2ID]
      DUnodes<-DUnodes[!DUnodes %in% common]
      DUnodes<-DUnodes[!duplicated(DUnodes)]
      DUgraph<-subKEGGgraph(DUnodes,bigGraph2)
      DUedges<-getKEGGedgeData(DUgraph)
      if(length(DUedges)>0){
        for(k in 1:length(DUedges)){
          if(strsplit(names(DUedges[k]),"~")[[1]][2] %in% D1nodes
             & strsplit(names(DUedges[k]),"~")[[1]][1] %in% U2ID){
            if(DUedges[[k]]@subtype[[1]]@name=="activation"
               | DUedges[[k]]@subtype[[1]]@name=="expression") DUn=DUn+1
            if(DUedges[[k]]@subtype[[1]]@name=="inhibition") DUp=DUp+1
          }
        }
      }
      ####DOWN+DOWN
      DID<-c(D1ID,D2ID)
      Dp<-0;Dn<-0;
      sharedPositive2<-Dshared[Dshared %in% bigGraphnodes]
      Dnodes<-DID[DID %in% bigGraphnodes2]
      Dgraph<-subKEGGgraph(Dnodes,bigGraph2)
      Dedges<-getKEGGedgeData(Dgraph)
      if(length(Dedges)>0){
        for(k in 1:length(Dedges)){
          if(strsplit(names(Dedges[k]),"~")[[1]][2] %in% D1nodes){
            if(Dedges[[k]]@subtype[[1]]@name=="activation"
               | Dedges[[k]]@subtype[[1]]@name=="expression") Dp=Dp+1
            if(Dedges[[k]]@subtype[[1]]@name=="inhibition") Dn=Dn+1
          }
        }
      }
      
      Dpositive=Dp+DUp-D1p
      Dnegtive=Dn+DUn-D1n
      
      positive<-Upositive+Dpositive
      negtive<-Unegtive+Dnegtive
      
      #### Part 2, estimate drug induced genomic residual effect
      ratio<-(length(sharedPositive)+length(sharedPositive2)+positive-negtive)
        /(length(U2ID)+length(D2ID)+length(Ushared)+length(Dshared))
      effect<-1-(1-ratio*doseRes2[j])*(1-(1-ratio)*doseRes1[j])*(1-doseRes1[i])
      effect<-effect-(1-(1-doseRes1[i])*(1-doseRes1[j]))
      
      resIC20<-rbind(resIC20,c(ratio,length(Dshared),length(sharedPositive), 
        length(sharedPositive2),positive,negtive,effect))
      rownames(resIC20)[length(resIC20[,1])]<-c(paste(names3[i],names3[j],sep="&"))
    }
  }
  
  ### the order of treatment is hypothesized to influence the results. To estimate the drug combination effect, we take the average of the estimations from 'B followed by A' cases  and 'A followed by B' cases.
  colnames(resIC20)<-c("ratio","shared Down","Upshared positive","downshared positive", 
    "positive effect(2to1)","negtive effect(2to1)","score(2to1)")
  write.csv(resIC20,file=resPath1)
  tempname<-strsplit(names(resIC20[,1]),"&")
  resIC20_2<-list()
  for(i in 1:length(tempname)){
    for(j in i:length(tempname)){
      if(i==j) next
      if(tempname[[i]][1]==tempname[[j]][2] 
         & tempname[[i]][2]==tempname[[j]][1]){
        resIC20_2<-rbind(resIC20_2,c(resIC20[i,1],resIC20[i,2],
                                     as.numeric(resIC20[i,7])+as.numeric(resIC20[j,7])))
        rownames(resIC20_2)[length(resIC20_2[,1])]<-rownames(resIC20)[i]         
      }
    }
  }
  colnames(resIC20_2)<-c("shared UP","shared Down","scores")
  write.csv(resIC20_2,file=resPath2)
}

