rm(list = ls())
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")


# 4. GSE103322 HNC
dataSet <- "GSE103322"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE103322_tg_data_list.RData")
tg_data_list <- GSE103322_tg_data_list 
Cell_Prop <- Cell_Prop_GSE103322

# set parameter
d=1
if(list_flag == T)	bulk <- tg_data_list[[d]][[1]]
if(list_flag == F)	bulk <- tg_data_list[[d]]
if(list_flag == T) true_p <- tg_data_list[[d]][[2]]
colnames(bulk) <- paste("ss",c(1:ncol(bulk)), sep="" )
colnames(true_p) <- paste("ss", c(1:ncol(true_p)), sep="")

#prepare high correlation data
ccc <- cor(t(bulk), t(true_p))
dim(ccc)
colnames(ccc)
B_mark <- names(ccc[order(-ccc[,3])[1:150],3])
T_mark <- names(ccc[order(-ccc[,2])[1:150],2])
Macro_mark <- names(ccc[order(-ccc[,8])[1:150],8])
Fibro_mark <- names(ccc[order(-ccc[,5])[1:150],5])
Epith_mark <- names(ccc[order(-ccc[,7])[1:150],7])
Mast_mark <- names(ccc[order(-ccc[,9])[1:150],9])

mark6 <- c(B_mark, T_mark, Macro_mark, Fibro_mark, Epith_mark, Mast_mark)

XXX <- bulk[mark6, ]
C <- matrix(0, nrow(XXX), 6)
rownames(C) <- rownames(XXX)
colnames(C) <- c("B","T","Macro","Fibro","EPith","Mast")
for(i in 1:ncol(C))
{	
	left <- 150*(i-1)+1
	right <- 150*i
	C[left:right, i] <- rep(1, 150)

}

bulk = XXX
ss_signature = C

#-------------------------------------------------run epic------------------------------
library(EPIC)


load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

EPIC_rmse <- list()
cancer_str <- dataSet

# just use gene which can produce protein
bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

# remove whole zero row
bulk <- rm_zero_row(bulk)

commonGene <- intersect(rownames(bulk),rownames(EPIC::TRef$refProfiles))
length(commonGene)
length(EPIC::TRef$sigGenes)
length(intersect(commonGene,EPIC::TRef$sigGenes))
commonSiga <- intersect(commonGene, EPIC::TRef$sigGenes)
commonSiga <- sort(commonSiga)
length(commonSiga)

sigGeneEpic <- EPIC::TRef$sigGenes
# the order in new_data_est is mess
sigGeneEpic <- sort(sigGeneEpic)

bulk1 <- bulk[commonSiga,]

signature <- EPIC::TRef$refProfiles[commonSiga,]
cell_marker <- extract_marker(signature,method="EPIC")
ss_signature <- unique_signature(signature, 0.05)
ss_signature <- rm_zero_row(ss_signature)
bulk1 <- bulk1[rownames(ss_signature), ]

out1 <- EPIC(bulk1)
print("finish EPIC")
Prop_EPIC <- out1$cellFraction[,1:7]

# do regression in order to get SS, then get bulk_est
SS <- do_regression_unique(Prop_EPIC, bulk1, ss_signature)
bulk_est <- SS %*% t(Prop_EPIC)

#R2
rmse_gene <- R2_two_mat(bulk_est, bulk1)

cell_rmse <- list()
for(i in 1:length(cell_marker)){
	cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
}
names(cell_rmse) <- names(cell_marker)

flag_plot <- T
if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_epic_cell_R2.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "R2", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()
}

EPIC_rmse[[1]] <- cell_rmse
names(EPIC_rmse)[1] <- cancer_str

save(EPIC_rmse, file = "EPIC_103322_list.RData")

#----------------------------ICTD------------------------------------------

library(EPIC)
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

ICTD_rmse <- list()
cancer_str <- dataSet

# just use gene which can produce protein
bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

# remove whole zero row
bulk <- rm_zero_row(bulk)

commonGene <- intersect(rownames(bulk),rownames(EPIC::TRef$refProfiles))
length(commonGene)
length(EPIC::TRef$sigGenes)
length(intersect(commonGene,EPIC::TRef$sigGenes))
commonSiga <- intersect(commonGene, EPIC::TRef$sigGenes)
commonSiga <- sort(commonSiga)
length(commonSiga)

sigGeneEpic <- EPIC::TRef$sigGenes
# the order in new_data_est is mess
sigGeneEpic <- sort(sigGeneEpic)

signature <- EPIC::TRef$refProfiles[commonSiga,]
cell_marker <- extract_marker(signature,method="EPIC")
ss_signature <- unique_signature(signature, 0.05)

ss_signature <- rm_zero_row(ss_signature)

ss_gene <- rownames(ss_signature)

############### NMF
data_t = log2(bulk + 1)			#NMF input 1
NMF_indi_all = ss_signature		#NMF input 2

X1=data_t[match(rownames(NMF_indi_all),rownames(data_t)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
###########

###########Parameter settings
theta=5 ##penalty parameter for constraints on NMF_indi_all
indiS_method="nonprdescent" ##the updating scheme for the structural constraints
iter=2000
alpha=beta=gamma=roh=0
roh=gamma=0.0001
UM=VM=NULL
qq=1
epslog=6
nPerm=2
initial_U=initial_V=NULL
mscale=1
###########


library(NMF)
set.seed(123456)
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/NMF/nmf.library.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/NMF/ini.R")
set.seed(123456)
###########Run the constrained qNMF
ttt1=qnmf_indisS_all_revise(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale)
#names(ttt1)
###########
U=ttt1$U[,c(1,2,4,3)]
V=ttt1$V[,c(1,2,4,3)]
U=ttt1$U
V=ttt1$V

SS <- do_regression_unique(V, X1, ss_signature)
bulk_est <- SS %*% t(V)

rmse_gene <- R2_two_mat(bulk_est, X1)

#extract marker corresponding rmse
cell_rmse <- list()
for(i in 1:length(cell_marker)){
	cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
}
names(cell_rmse) <- names(cell_marker)


flag_plot <- T
if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_ICAD_cell_R2.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "R2", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()
}

ICTD_rmse[[1]] <- cell_rmse
names(ICTD_rmse)[1] <- cancer_str

save(ICTD_rmse, file = "ICTD_103322_list.RData")


#--------------------------ICTD-self-marker-----------------------------------


load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

ICTD_rmse <- list()
cancer_str <- dataSet

# just use gene which can produce protein
bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

# remove whole zero row
bulk <- rm_zero_row(bulk)
ss_signature <- ss_signature[rownames(bulk),]
cell_marker <- extract_marker_ICTD(ss_signature)


############### NMF
data_t = log2(bulk + 1)			#NMF input 1
NMF_indi_all = ss_signature		#NMF input 2

X1=data_t[match(rownames(NMF_indi_all),rownames(data_t)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
###########

###########Parameter settings
theta=5 ##penalty parameter for constraints on NMF_indi_all
indiS_method="nonprdescent" ##the updating scheme for the structural constraints
iter=2000
alpha=beta=gamma=roh=0
roh=gamma=0.0001
UM=VM=NULL
qq=1
epslog=6
nPerm=2
initial_U=initial_V=NULL
mscale=1
###########


library(NMF)
set.seed(123456)
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/NMF/nmf.library.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/NMF/ini.R")
set.seed(123456)
###########Run the constrained qNMF
ttt1=qnmf_indisS_all_revise(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale)
#names(ttt1)
###########
U=ttt1$U[,c(1,2,4,3)]
V=ttt1$V[,c(1,2,4,3)]
U=ttt1$U
V=ttt1$V

#1
SS <- do_regression_unique(V, X1, ss_signature)
bulk_est <- SS %*% t(V)
rmse_gene <- R2_two_mat(bulk_est, X1)

#2
#SS <- correspond_matrix(U, ss_signature)
#bulk_est <- SS %*% t(V)
######bulk_est <- U %*% t(V)
#rmse_gene <- R2_two_mat_withIntercept(bulk_est, X1)


#extract marker corresponding rmse
cell_rmse <- list()
for(i in 1:length(cell_marker)){
	cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
}
names(cell_rmse) <- names(cell_marker)


flag_plot <- T
if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_ICAD_cell_R2_directRegression.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "R2", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()
}

ICTD_rmse[[1]] <- cell_rmse
names(ICTD_rmse)[1] <- cancer_str

save(ICTD_rmse, file = "ICTD_103322_list.RData")

































##########################################################
 #
 #    function
 #########################################################

filter_gene_name <- function(data_t){
# 1. extra special row
loca <- c()
for(i in 1:nrow(data_t)){
	tmp <- unlist(strsplit(rownames(data_t)[i], '|', fixed = T))
	if(length(tmp) > 1) loca <- c(loca, i)
}
data_sub <- data_t[loca,]

# 2. change name
gene <- c()
for(i in 1:nrow(data_sub)){
	tmp <- unlist(strsplit(rownames(data_sub)[i], '|', fixed = T))
	gene <- c(gene, tmp)
}
gene_2col <- matrix(gene, 2, length(gene)/2)
gene_name <- gene_2col[2,]
rownames(data_sub) <- gene_name

data_ttt <- data_sub

return(data_ttt)

}


getFractions.Abbas <- function(XX,YY,w=NA){
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0,ncol(XX)))
      tmp.XX=tmp.XX[,-ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0,ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX,YY,intercept=F) else tmp=lsfit(tmp.XX,YY,w,intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0,ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}


TIMER <- function(cc = cancer_str){
	if(cc=='skcm')cc.type='06A' else cc.type='01A'

	##----- setup parameters and establish the output file -----##
	signature.genes=c('CD19','TRAT1','CD8B','CCR3','CD163','CCL17')
	names(signature.genes)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

	##----- load and process gene expression data -----##
	file_str <- paste("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA-", toupper(cancer_str), "_FPKM_T.RData", sep="")
	data_t <- get(load(file_str))
	data_ttt <- filter_gene_name(data_t)
	dd <- data_ttt
	#dd <- get(load(file_str) )
	dd=as.matrix(dd)
	mode(dd)='numeric'	
	#if(!cc %in% c('gbm','ov','esca','stad'))dd=dd*1e6   ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
	tmp=strsplit(rownames(dd),'\\|')
	tmp=sapply(tmp,function(x)x[[1]])
	tmp.vv=which(nchar(tmp)>1)
	rownames(dd)=tmp
	dd=dd[tmp.vv,]

	##----- load immune marker genes from Abbas et al., 2005 -----##
	tmp=read.csv("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/IRIS-marker-gene.txt", ,header=T,sep='\t',stringsAsFactors=F) # * add / before data
	marker.list=tmp[,1]
	names(marker.list)=tmp[,7]
	names(marker.list)=gsub(' ','_',tmp[,7])
	names(marker.list)=gsub('Dendritic_Cell','DC',names(marker.list))
	names(marker.list)=gsub('Neutrophil','Neutrophils',names(marker.list))
	gsub('Cell','cell',names(marker.list))->names(marker.list)

	##----- load reference data of sorted immune cells -----##
	load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/HPCTimmune.Rdata") # * add / before data

	##----- load and process tumor purity data -----##
	AGP=read.table(paste('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/AGP/AGP-',cc,'.txt',sep=''),sep='\t',header=T) # * add / before datas
	AGP=AGP[which(AGP[,'PoP']>0.01),]
	tmp=strsplit(rownames(HPCT.immune),';')
	AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
	names(AffyIDtoGenes)=sapply(tmp,function(x)x[[2]])
	marker.list.genes=AffyIDtoGenes[marker.list]

	##----- function to edit TCGA ID, with the option of keeping the first num.res fields -----##
	getID <- function(sID,num.res=3){
 		mm=c()
  		for(id in sID){
   			tmp=unlist(strsplit(id,'-'))
    		if(length(tmp)==1){
     	 		tmp=unlist(strsplit(id,'\\.'))
    		}
    		ll='TCGA'
    		for(j in 2:num.res){
     			ll=paste(ll,tmp[j],sep='-')
   		 	}
    		mm=c(mm,ll)
 		}
  		return(mm)
	}
	rownames(AGP)=getID(AGP[,1],4)
	colnames(dd)=getID(colnames(dd),4)


	##----- Select single reference samples of pre-selected immune cell types -----##
	B_cell=362:385
	T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
	T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
	NK=328:331
	Neutrophil=344:361
	Macrophage=66:80
	DC=151:238
	curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

	curated.cell.types=colnames(curated.ref)
	names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))

	load('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/curated.ref.genes.Rdata')

	##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
	RemoveBatchEffect <- function(){
  		library(sva)
  		tmp.dd=as.matrix(dd)
  		tmp=sapply(strsplit(rownames(dd),'\\|'),function(x)x[[1]])
  		rownames(tmp.dd)=tmp
  		tmp.dd=tmp.dd[which(nchar(tmp)>1),]
  		tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  		N1=ncol(tmp.dd)
  		tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,])
  		tmp.dd=as.matrix(tmp.dd)
  		mode(tmp.dd)='numeric'
  		N2=ncol(curated.ref.genes)
  		tmp.batch=c(rep(1,N1),rep(2,N2))
  		tmp.dd0=ComBat(tmp.dd,tmp.batch,c())
  		dd.br=tmp.dd0[,1:N1]
  		curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)]
  		tmp0=c()
  		for(kk in unique(names(curated.cell.types))){
   		 tmp.vv=which(names(curated.cell.types)==kk)
   		 tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  		}
  		curated.ref.genes.agg.br=tmp0
  		colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  		#rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  		return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
	}

	tmp=RemoveBatchEffect()
	dd.br=tmp$dd
	curated.ref.genes.br=tmp$rr
	curated.ref.genes.agg.br=tmp$rrg


	##----- function to calculate the residuals from regression -----##
	fn <- function(beta0,XX,Y)return(log(sum(abs(Y-XX%*%beta0))))

	##----- function to select genes with expression values negatively correlated with tumor purity -----##
	getPurityGenes <- function(dd,AGP,thr.p=0.05,thr.c=0,mode='env'){
 	 tmp.ss=intersect(colnames(dd),rownames(AGP))
  	if(length(tmp.ss)==0){
    	colnames(dd)=getID(colnames(dd))
   		tmp.ss=intersect(colnames(dd),rownames(AGP))
  	}
  	tmp.dd=dd[,tmp.ss]	#original
 	#tmp.dd <- dd    #w
  	tmp=lapply(rownames(tmp.dd),function(x)cor.test(tmp.dd[x,],as.numeric(AGP[colnames(tmp.dd),2]),method='s'))
  	tmp.pp=sapply(tmp,function(x)x$p.value)
  	tmp.cor=sapply(tmp,function(x)x$estimate)
  	names(tmp.pp)=names(tmp.cor)=rownames(dd)
  	if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
 	if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
 	return(vv)
	}

	##----- selection genes negatively correlated with purity and overlap with immune marker genes -----##
	vv.t=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= -0.2)
	vv.t=intersect(vv.t,rownames(curated.ref.genes.agg.br))
	vv=intersect(vv.t,marker.list.genes)

	##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
	RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
 	 ## removes upper thr.q quantile for every reference feature
 	 remove.vv=c()
 	 for(i in 1:ncol(ref.dd)){
  	  tmp=quantile(ref.dd[vv,i],thr.q)[1]
  	  tmp.vv=which(ref.dd[vv,i]>tmp)
 	  remove.vv=c(remove.vv,tmp.vv)
	 }
 	 remove.vv=unique(remove.vv)
 	 return(vv[-remove.vv])
	}

	##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
	tmp.diff=sum(sum(abs(cor(curated.ref.genes.agg.br[vv,],method='p')-cor(curated.ref.genes.agg.br[vv,],method='s'))))

	if(tmp.diff>= -10000){
 		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4])
 		 vv=vv0
	}

	cat("Number of genes inversely correlated with purity is ",length(vv.t),'\n\n',sep='',file='output-statistics.txt',append=T)
	cat("Number of immune genes inversely correlated with purity is ",length(vv),'\n\n',sep='',file='output-statistics.txt',append=T)

	##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
	tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
	n.immune=length(intersect(marker.list.genes,tmp.ss0))
	cat("Test if immune genes are enriched for inverse correlation with purity: \n\n",file='output-statistics.txt',append=T)
	sink(file='output-statistics.txt',append=T);print(fisher.test(matrix(c(length(vv),length(vv.t)-length(vv),n.immune,length(tmp.ss0)-n.immune),2,2)));sink()

	##----- function to process deconvolution method in batch -----##
	BatchFractions <- function(XX,YYd){
  		Fmat=c()
  		for(i in 1:ncol(YYd)){
  			YY=YYd[,i]
   			tmp.F=getFractions.Abbas(XX,YY)
    		#tmp.F=getFractions.Optim(XX,YY)
    		Fmat=rbind(Fmat,tmp.F)
  		}
  		rownames(Fmat)=colnames(YYd)
  		colnames(Fmat)=colnames(XX)
  		return(Fmat)
	}

	##----- perform batch deconvolution -----##
	XX=curated.ref.genes.agg.br[vv,c(-4)]	#chang note: delete NK cell during deconvolution
	YYd=dd.br[vv,]
	Fmat=BatchFractions(XX,YYd)


	##----- CD4 and CD8 T cells are likely to be similar, resulting in colinearity. Codes below are procedures to remove outlier genes that may result in colinearity until the two covariates are linearly separable. -----##
	if(cor(Fmat[,2],Fmat[,3])<= -0.2){
  		if(tmp.diff>=1){
    		tmp.cor=c()
    		thr.qlist=c(0.99)
    		for(tq in thr.qlist){
     		 vv=intersect(vv.t,marker.list.genes)
    		 vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4],tq)
    		 vv=vv0
     		 XX=curated.ref.genes.agg.br[vv,]
     		 YYd=dd.br[vv,]
     		 tmp.Fmat=BatchFractions(XX,YYd)
     		 tmp.cor=c(tmp.cor,cor(tmp.Fmat[,2],tmp.Fmat[,3],method='s'))
   			}
   			tmp.vv=which.max(tmp.cor)
    		tq=thr.qlist[tmp.vv]
    		vv=intersect(vv.t,marker.list.genes)
    		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,c(-4)],tq)
    		vv=vv0
    		XX=curated.ref.genes.agg.br[vv,c(-4)]
    		YYd=dd.br[vv,]
    		Fmat=BatchFractions(XX,YYd)
    		Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
    		rownames(Fmat0.p)=getID(rownames(Fmat0.p))
  		}
	}

	while(cor(Fmat[,2],Fmat[,3])<=-0.3){
 	 if(length(vv)<=50)break
	 vv=vv[-as.numeric(names(table(apply(dd[vv,],2,which.max))))]
  	 XX=curated.ref.genes.agg.br[vv,c(-4)]
     #XX=XX[,!colnames(XX) %in% tmp.remove]
     YYd=dd.br[vv,]
     Fmat=BatchFractions(XX,YYd)
	}

	Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
	rownames(Fmat0.p)=getID(rownames(Fmat0.p))

	TIMER_fraction <- Fmat
	X <- YYd
	S <- XX #negative value
	return(list(Timer_fraction = TIMER_fraction, X = X, S = S))

}




RMSE_one_vector <- function(a){
	rmse <- sqrt(sum(a^2)/length(a))
	return(rmse)
}

RMSE_two_vector <- function(a, b){
	diff <- a - b
	rmse <- sqrt(sum(diff^2)/length(diff))
	return(rmse)
}

R2_two_vector <- function(predict, actual){
	#R2 <- 1 - sum( (actual-predict )^2 ) / sum( (actual-mean(actual) )^2  ) 
	R2 <- 1 - sum( (actual-predict )^2 ) / sum( actual^2 ) 
	return(R2)
}

R2_two_vector_withIntercept <- function(predict, actual){
	R2 <- 1 - sum( (actual-predict )^2 ) / sum( (actual-mean(actual) )^2  ) 
	return(R2)
}

RMSE_one_mat <- function(aaa){
	
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		vc[i] <- RMSE_one_vector(tmp_a)

	}
	rownames(vc) <- rownames(aaa)
	
	return(vc)
}

RMSE_two_mat <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- RMSE_two_vector(tmp_a, tmp_b)

	}
	rownames(vc) <- rownames(aaa)
	
	return(vc)
}

R2_two_mat <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- R2_two_vector(tmp_a, tmp_b)
	}
	rownames(vc) <- rownames(aaa)
	return(vc)
}

R2_two_mat_withIntercept <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- R2_two_vector_withIntercept(tmp_a, tmp_b)
	}
	rownames(vc) <- rownames(aaa)
	return(vc)
}



 do_regression_unique <- function(proportion=Prop_EPIC[,1:7], data=bulk1, indicate=ss_signature){
 	print("unique regression")
 	if(nrow(data) != nrow(indicate)) stop("number of gene in indicate NOT match!")
 	if(ncol(proportion) != ncol(indicate)) stop("number of sample in indicate NOT match!")
 	P_2nd <- matrix(0, nrow(indicate), ncol(indicate))
 	rownames(P_2nd) <- rownames(indicate)
	n_gene <- nrow(data)
	n_sample <- ncol(data)	
	for(i in 1:n_gene){
		#print(i)
		cell <- which(indicate[i, ] == 1)
		#if(length(cell) > 1) stop("return 2 cell type!")
	 	pp <- proportion[, cell]
		y <- data[i, ]
		ddd <- cbind(pp, y)
		ddd <- as.data.frame(ddd)
		coeff <- lm(y~pp + 0, ddd)	#NOTE : the order is matter!
		P_2nd[i, cell] <- coeff[[1]]
	}	
	return(P_2nd)
 }

extract_marker <- function(aaa, method){
	# mmmm is a list to store cell marker, the name of each element is cell type,
	# 
	cutoff <- 0.05
	mmm <- matrix(0, nrow(aaa), ncol(aaa))
	rownames(mmm) <- rownames(aaa)
	colnames(mmm) <- colnames(aaa)
	for(i in 1:nrow(aaa)){
		vv <- matrix(NA, 1, ncol(aaa))
		vv <- aaa[i, ]

		gene_sd <- sd(vv) * sqrt( (length(vv)-1) / (length(vv))  )
		gene_mean <- mean(vv)

		#z-score calculation
		z <- (vv - gene_mean) / gene_sd

		p_yellow <- pnorm(z)
		p_blue <- 1 - p_yellow
		#which(p_blue < cutoff)
		mmm[i, which(p_blue < cutoff)] <- 1
	}

	#find zero rows
	#mmm <- extract_marker(aaa)
	#row_sub <- apply(mmm, 1, function(row) all( row == 0)) #return logistic value(T/F)
	#Zero <- which(row_sub == T)
	#length(Zero) 
	#length(which(mmm[, 1] != 0))
	#length(which(mmm[, 2] != 0) )
	#length(which(mmm[, 3] != 0) )
	#length(which(mmm[, 4] != 0) )
	#length(which(mmm[, 5] != 0) )
	#length(which(mmm[, 6] != 0) )

	if(method == "TIMER"){
		B_mark <- rownames(mmm)[which(mmm[, 1] != 0) ]
		CD4_mark <- rownames(mmm)[which(mmm[, 2] != 0) ]
		CD8_mark <-rownames(mmm)[which(mmm[, 3] != 0) ]
		Netro_mark <- rownames(mmm)[which(mmm[, 4] != 0) ]
		Macro_mark <- rownames(mmm)[which(mmm[, 5] != 0) ]
		Dc_mark <- rownames(mmm)[which(mmm[, 6] != 0) ]
	}
	#epic
	if(method == "EPIC"){
		B_mark <- rownames(mmm)[which(mmm[, 1] != 0) ]
		CAFs_mark <- rownames(mmm)[which(mmm[, 2] != 0) ]
		CD4_mark <-rownames(mmm)[which(mmm[, 3] != 0) ]
		CD8_mark <- rownames(mmm)[which(mmm[, 4] != 0) ]
		endoth_mark <- rownames(mmm)[which(mmm[, 5] != 0) ]
		Macro_mark <- rownames(mmm)[which(mmm[, 6] != 0) ]
		NK_mark <- rownames(mmm)[which(mmm[, 7] != 0) ]
	}
	#ciber
	if(method == "CIBER"){
		specific_marker <- list()
		for(i in 1:ncol(mmm)){
			specific_marker[[i]] <- rownames(mmm)[which(mmm[, i] != 0) ]
			names(specific_marker)[i] <- colnames(mmm)[i]
		}
	}


	#return
	if(method == "TIMER"){
		return(list(B_mark = B_mark, CD4_mark = CD4_mark, CD8_mark = CD8_mark, Netro_mark = Netro_mark, Macro_mark =Macro_mark, DC_mark = Dc_mark))
	}
	if(method == "EPIC"){
		return(list(B_mark = B_mark, CAFs_mark = CAFs_mark, CD4_mark = CD4_mark, CD8_mark = CD8_mark, endoth_mark =endoth_mark, Macro_mark = Macro_mark, NK_mark=NK_mark))
	}
	if(method == "CIBER"){
		return(specific_marker)
	}


}

extract_marker_ICTD <- function(C){
	specific_marker <- list()
	for(i in 1:ncol(C)){
		specific_marker[[i]] <- rownames(C)[which(C[, i] != 0)]
		names(specific_marker)[i] <- colnames(C)[i]
	}

	return(specific_marker)
}

unique_signature <- function(aaa, cutoff=0.05){

	#cutoff <- 0.05
	mmm <- matrix(0, nrow(aaa), ncol(aaa))
	rownames(mmm) <- rownames(aaa)
	colnames(mmm) <- colnames(aaa)
	for(i in 1:nrow(aaa)){
		vv <- matrix(NA, 1, ncol(aaa))
		vv <- aaa[i, ]

		gene_sd <- sd(vv) * sqrt( (length(vv)-1) / (length(vv))  )
		gene_mean <- mean(vv)

		#z-score calculation
		z <- (vv - gene_mean) / gene_sd

		p_yellow <- pnorm(z)
		p_blue <- 1 - p_yellow
		#which(p_blue < cutoff)
		mmm[i, which(p_blue < cutoff)] <- 1
	}

	return(mmm)
}


rm_zero_row <- function(bulk){

	row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) bulk <- bulk
	if(length(Zero) > 0)  bulk <- bulk[-Zero,]

	return(bulk)
}

correspond_matrix <- function(sss=U, ccc=ss_signature){
	mmm <- matrix(0, nrow(sss), ncol(sss))
	for(i in 1:nrow(ccc)){
		cell <- which(ccc[i,] == 1)
		mmm[i, cell] <- sss[i, cell]
	}
	rownames(mmm) <- rownames(sss)
	colnames(mmm) <- colnames(sss)

	return(mmm)
}






