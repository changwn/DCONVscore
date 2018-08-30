#
#
#
library(EPIC)
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

#---------------------------------------------------------------------------
# load data, which is ttt
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/RNAseq/coadRNAseq.RData")

# remove same gene
ttt1 <- ttt[unique(rownames(ttt)),]

bulk <- ttt1
# out <- EPIC(bulk, scaleExprs=F)
out <- EPIC(bulk)
names(out)

commonGene <- intersect(rownames(bulk),rownames(EPIC::TRef$refProfiles))
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
bulk1 <- new_data
# out1 <- EPIC(bulk1, scaleExprs=F)
out1 <- EPIC(bulk1)
Prop_EPIC <- out1$cellFraction[,1:7]

# names(out1)
# If scaleExprs is false, we find the results of cellFraction of out and out1 is exactly same

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
#rownames(ccc) <- colnames(ccc)
dd <- diag(ccc)
dd_copy <- dd
Dscore_whole <- mean(dd)

cal_Dscore <- function(proportion=prop_base, data=bulk1){
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
	dd_copy <- dd
	Dscore <- mean(dd)
	return(Dscore)
}

# choice 1:Find top10 genes which has high correlation with the original data 
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
Dscore_init <- cal_Dscore(prop_base, bulk1)
}

# choice 2: sort the gene based on correlation and then extract top 10 gene
dd_copy_order <- sort(dd_copy, decreasing = T)
base <- names(dd_copy_order[1:10])
out_base <- EPIC(bulk1, sigGenes = base)
prop_base <- out_base$cellFraction[, 1:7]
Dscore_init <- cal_Dscore(prop_base, bulk1)


#remainSet <- setdiff(sigGeneEpic, base)
remainSet <- names(dd_copy_order[11:98])
top_gene_add <- base
Dscore_Yaxis <- c()
increase_gene <- c()
for(i in 1:length(remainSet)){
#for(i in 1 : 5){
	top_gene_add <- union(top_gene_add, remainSet[i])
	out_ep <- EPIC(bulk1, sigGenes=top_gene_add)
	prop_add <- out_ep$cellFraction[,1:7]
	Dscore_add <- cal_Dscore(prop_add, bulk1)
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
# results : the order put gene in the "top_gene_add" will change the final gene list and score.
# should check the rank within the interation

x <- c(11:98)
wholeName <- matrix(NA,1,88)
wholeName[match(increase_gene, remainSet)] <- increase_gene
plot(x, Dscore_Yaxis, type = 'l', main = "coad, score")
text(x, Dscore_Yaxis, wholeName, cex=0.8,  col = "red", srt = 30)
gene_char <- paste(setdiff(top_gene_add, base), sep=" ", collapse=",")
text(55, 0.26, gene_char, cex=0.6,  col = "blue")
base_char <- paste(base, sep="", collapse=",")
text(38, 0.265, base_char, cex=0.6, col="black")

# --------------------------------------------------------------------
# evaluation on single cell simulated data
# melanoma data
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE72056_tg_data_list.RData")

# using 28 top_gene_add
storage1_top <- list()
names(storage1) <- names(Cell_Prop_GSE72056)
for(i in 1:length(Cell_Prop_GSE72056)){

	bulk1 <- GSE72056_tg_data_list[[i]][[1]]	
	#bulk1 <- bulk1[top_gene_add,]

	out1 <- EPIC(bulk1, sigGenes = top_gene_add)
	prop_sc <- out1$cellFraction[, 1:7]
	prop_true <- Cell_Prop_GSE72056[[i]]
	prop_true <- t(prop_true)
	corr_sc <- cor(prop_sc, prop_true)
	corr_sc_top <- corr_sc
	storage1_top[[length(storage1_top)+1]] <- corr_sc_top
}



# using all sigGeneEpic
storage1_all <- list()
names(storage1_all) <- names(Cell_Prop_GSE72056)
for(i in 1:length(Cell_Prop_GSE72056)){

	bulk1 <- GSE72056_tg_data_list[[i]][[1]]	
	#bulk1 <- bulk1[top_gene_add,]

	out1 <- EPIC(bulk1)
	prop_sc <- out1$cellFraction[, 1:7]
	prop_true <- Cell_Prop_GSE72056[[i]]
	prop_true <- t(prop_true)
	corr_sc <- cor(prop_sc, prop_true)
	corr_sc_all <- corr_sc
	storage1_all[[length(storage1_all)+1]] <- corr_sc_all
}

diff <- corr_sc_top - corr_sc_all
# results show few gene signature still produce high correlation,
# excepting T cell(CD4, CD8 T cell)



