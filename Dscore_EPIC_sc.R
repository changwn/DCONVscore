
library(EPIC)
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")


# 1. GSE72056 melanoma
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE72056_tg_data_list.RData")
tg_data_list <- GSE72056_tg_data_list 
Cell_Prop <- Cell_Prop_GSE72056
# 2. GSE70630 oligodendroglioma
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE70630_tg_data_list.RData")
tg_data_list <- GSE70630_tg_data_list 
Cell_Prop <- Cell_Prop_GSE70630
# 3. GSE89567 astrocytoma
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE89567_tg_data_list.RData")
tg_data_list <- GSE89567_tg_data_list 
Cell_Prop <- Cell_Prop_GSE89567
# 4. GSE103322 HNC
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE103322_tg_data_list.RData")
tg_data_list <- GSE103322_tg_data_list 
Cell_Prop <- Cell_Prop_GSE103322


# 5. GSE75688 breast cancer
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_29_singleCellSimulated/GSE75688_tg_data_list.RData")
tg_data_list <- GSE75688_tg_data_list 
Cell_Prop <- Cell_Prop_GSE75688
# 6. GSE81861 colorectal tumors
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_29_singleCellSimulated/GSE81861_tg_data_list.RData")
tg_data_list <- GSE81861_tg_data_list 
Cell_Prop <- Cell_Prop_GSE81861


i=1
bulk <- tg_data_list[[i]][[1]]
#bulk <- tg_data_list[[i]]
out <- EPIC(bulk)

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
S <- EPIC::TRef$refProfiles[sigGeneEpic,]

new_data <- bulk[commonSiga,]
# delete 0 row
row_sub <- apply(new_data, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk1 <- new_data
if(length(Zero) > 0)  bulk1 <- new_data[-Zero,]
# out1 <- EPIC(bulk1, scaleExprs=F)
out1 <- EPIC(bulk1)
Prop_EPIC <- out1$cellFraction[,1:7]

# Use predicted proportion(P matrix) to find how much it can explain from data(X matrix),
##----- Constrained regression method implemented in Abbas et al., 2009 -----##
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


P_2nd <- c()
n_gene <- nrow(bulk1)
n_sample <- ncol(bulk1)
for(i in 1:n_gene){
	coeff <- getFractions.Abbas(Prop_EPIC, t(bulk1)[,i])	# first,use complete proportion to get rank of gene
	g_tmp <- rownames(bulk1)[i]
	P_2nd <- rbind(P_2nd, coeff)
	rownames(P_2nd)[i] <- g_tmp 
}
bulk_est <- P_2nd %*% t(Prop_EPIC)
ccc <- cor(t(bulk_est), t(bulk1))	# want to get rank of gene-wise
rownames(ccc) <- colnames(ccc)
dd <- diag(ccc)
# absolute value
dd <- abs(dd)
dd_copy <- dd
Dscore_whole <- mean(dd)

cal_Dscore_abs <- function(proportion=prop_base, data=bulk1){
	P_2nd <- c()
	n_gene <- nrow(data)
	n_sample <- ncol(data)
	for(i in 1:n_gene){
		coeff <- getFractions.Abbas(proportion, t(data)[, i])
		P_2nd <- rbind(P_2nd, coeff)
	}
	bulk_est <- P_2nd %*% t(proportion)
	# then, how to evaluation the similarity of two matrix (bulk and bulk_est)?
	# try sample-wise correlation 
	ccc <- cor(t(bulk1), t(bulk_est))	#??? still gene-wise
	dd <- diag(ccc)
	dd <- abs(dd)
	dd_copy <- dd
	Dscore <- mean(dd)
	return(Dscore)
}

# choice 1
while(F){
topN <- 10
top_gene <- c()
for(i in 1:topN){
	gene_tmp <- names(dd[which(dd_copy == max(dd_copy))])
	assign(paste("top", i, sep=""), gene_tmp)
	dd_copy[gene_tmp] <- 0
	top_gene <- c(top_gene, gene_tmp)
}
base <- top_gene
out_base <- EPIC(bulk1, sigGenes = base)
prop_base <- out_base$cellFraction[,1:7]
Dscore_init <- cal_Dscore_abs(prop_base, bulk1)
}

# choice 2: sort the gene based on correlation and then extract top 10 gene
dd_copy_order <- sort(dd_copy, decreasing = T)
base <- names(dd_copy_order[1:10])
out_base <- EPIC(bulk1, sigGenes = base)
prop_base <- out_base$cellFraction[, 1:7]
Dscore_init <- cal_Dscore_abs(prop_base, bulk1)

# score of explaination level of original data
while(F){
#remainSet <- setdiff(sigGeneEpic, top_gene)
#top_gene_add <- top_gene
remainSet <- names(dd_copy_order[11:98])
top_gene_add <- base
Dscore_Yaxis <- c()
increase_gene <- c()
for(i in 1:length(remainSet)){
#for(i in 15:length(remainSet)){
	top_gene_add <- union(top_gene_add, remainSet[i])
	out_ep <- EPIC(bulk1, sigGenes=top_gene_add)
	prop_add <- out_ep$cellFraction[,1:7]
	Dscore_add <- cal_Dscore_abs(prop_add, bulk1)
	Dscore_Yaxis[i] <- Dscore_add
	print(Dscore_add)
	if(Dscore_add < Dscore_init){
		top_gene_add <- top_gene_add[1:length(top_gene_add)-1]
		print("less, delete")
	}else{
		Dscore_init <- Dscore_add
		print("large, change init value")
		increase_gene <- c(increase_gene, remainSet[i])
	}
}
print(top_gene_add)
}

x <- c(1:88)
wholeName <- matrix(NA,1,88)
wholeName[match(increase_gene, remainSet)] <- increase_gene
plot(x, Dscore_Yaxis, type = 'l', main = "GSE103322 HNC, score")
text(x, Dscore_Yaxis, wholeName, cex=0.8,  col = "red", srt = 30)
base_char <- paste(base, sep="", collapse=",")
text(35, 0.755, base_char, cex=0.7, col="black")
gene_char <- paste(setdiff(top_gene_add,base), sep=" ", collapse=",")
text(40, 0.751, gene_char, cex=0.7,  col = "blue")


# evaluation
# using 17 top_gene_add
storage1_top <- list()
for(i in 1:length(Cell_Prop)){

	#bulk1 <- tg_data_list[[i]][[1]]	
	bulk1 <- tg_data_list[[i]]

	out1 <- EPIC(bulk1, sigGenes = top_gene_add)
	#prop_sc <- out1$cellFraction[, 1:7]
	prop_sc <- out1$cellFraction[, 1:8]
	prop_true <- Cell_Prop[[i]]
	prop_true <- t(prop_true)
	corr_sc <- cor(prop_sc, prop_true)
	corr_sc_top <- corr_sc
	storage1_top[[length(storage1_top)+1]] <- corr_sc_top
}
names(storage1_top) <- names(Cell_Prop)

# using all sigGeneEpic
storage1_all <- list()
for(i in 1:length(Cell_Prop)){

	#bulk1 <- tg_data_list[[i]][[1]]	
	bulk1 <- tg_data_list[[i]]

	out1 <- EPIC(bulk1)
	#prop_sc <- out1$cellFraction[, 1:7]
	prop_sc <- out1$cellFraction[, 1:8]
	prop_true <- Cell_Prop[[i]]
	prop_true <- t(prop_true)
	corr_sc <- cor(prop_sc, prop_true)
	corr_sc_all <- corr_sc
	storage1_all[[length(storage1_all)+1]] <- corr_sc_all
}
names(storage1_all) <- names(Cell_Prop)











