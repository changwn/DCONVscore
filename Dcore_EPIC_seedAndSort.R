
# use common genes of EPIC and CIBERSORT as seed (around 29 genes),
# sort all signature genes in CIBERSORT (547 genes) from high to low based on explaination idea.
# Finally, evaluate and plot the correlation score as gene list accumulated.

library(EPIC)
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")

# 1. GSE72056 melanoma
dataSet = "GSE72056"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE72056_tg_data_list.RData")
tg_data_list <- GSE72056_tg_data_list 
Cell_Prop <- Cell_Prop_GSE72056

# set parameter
d=1
if(list_flag == T)	bulk <- tg_data_list[[d]][[1]]
if(list_flag == F)	bulk <- tg_data_list[[d]]

sigGeneEpic <- EPIC::TRef$sigGenes
sigGeneEpic <- sort(sigGeneEpic)
sigGeneCiber <- rownames(SigMat_CIBER)
sigGeneCiber <- sort(sigGeneCiber)
length(sigGeneEpic)
length(sigGeneCiber)
length(intersect(sigGeneCiber, sigGeneEpic))

commonGene <- intersect(rownames(bulk), sigGeneCiber)
length(commonGene)

#length(intersect(commonGene,EPIC::TRef$sigGenes))
#commonSiga <- intersect(commonGene, EPIC::TRef$sigGenes)
#commonSiga <- sort(commonSiga)
#length(commonSiga)

#same_EandC <- intersect(sigGeneEpic, sigGeneCiber)
#length(same_EandC)

new_data <- bulk[commonGene,]
# delete 0 row
row_sub <- apply(new_data, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk1 <- new_data
if(length(Zero) > 0)  bulk1 <- new_data[-Zero,]
# out1 <- EPIC(bulk1, scaleExprs=F)
out1 <- EPIC(bulk1)
Prop_EPIC <- out1$cellFraction[,1:7]

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

# get TOP 10: sort the gene based on correlation and then extract top 10 gene
dd_copy_order <- sort(dd_copy, decreasing = T)
base <- names(dd_copy_order[1:10])
out_base <- EPIC(bulk1, sigGenes = base)
prop_base <- out_base$cellFraction[, 1:7]
Dscore_init <- cal_Dscore_abs(prop_base, bulk1)



# initial correlation score
out_base <- EPIC(bulk, sigGenes = same_EandC)
prop_base <- out_base$cellFraction[, 1:8]
prop_base <- sum_cd4_cd8(prop_base)
prop_true <- t(Cell_Prop[[d]])
corr_base <- cor(prop_base, prop_true)

if(dataSet == "GSE72056") Dscore_mat <- matrix(NA, length(remainSet), 8)
if(dataSet == "GSE70630") Dscore_mat <- matrix(NA, length(remainSet), 1)
if(dataSet == "GSE89567") Dscore_mat <- matrix(NA, length(remainSet), 1)
if(dataSet == "GSE103322") Dscore_mat <- matrix(NA, length(remainSet), 9)
if(dataSet == "GSE75688") Dscore_mat <- matrix(NA, length(remainSet), 7)
if(dataSet == "GSE81861") Dscore_mat <- matrix(NA, length(remainSet), 7)

Dscore_mat <- matrix(NA, 1, 8)


i=1
if(dataSet == "GSE72056"){	
	Dscore_mat[i, 1] <- corr_base[1, 1]	#B cell
	Dscore_mat[i, 2] <- corr_base[9, 2]	#T
	Dscore_mat[i, 3] <- corr_base[8, 3]	#tumor
	Dscore_mat[i, 4] <- corr_base[6, 4]	#macrophage
	Dscore_mat[i, 5] <- corr_base[2, 5]	#Fibroblast
	Dscore_mat[i, 6] <- corr_base[5, 6]	#endothelial
	Dscore_mat[i, 7] <- mean(Dscore_mat[i, 1:6])	#average
	tmp <- (sum(Dscore_mat[i,1:6]) - Dscore_mat[i, 3]) / 5
	Dscore_mat[i, 8] <- tmp  			#ave-5
	if(i==1) colnames(Dscore_mat) <- c("B","T","tumor","macrophage","Fibroblast","endothelial","ave_all","ave_5")

}

remainSet <- setdiff(sigGeneCiber, same_EandC)

# First, use explanation idea to give the order of all gene,
# then, use this order to calculate the correlation scoreby adding one gene each time.
# 1. if score not decend, means that above order not work?
# 2. the group of gene (row space) cannot be evaluate by adding one each time?
# 3. the cibersort list not work?

new_data <- bulk[sigGeneCiber,]
# delete 0 row
row_sub <- apply(new_data, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk1 <- new_data
if(length(Zero) > 0)  bulk1 <- new_data[-Zero,]
# out1 <- EPIC(bulk1, scaleExprs=F)
out1 <- EPIC(bulk1)
Prop_EPIC <- out1$cellFraction[,1:7]




#-----------------------------
#function

# add a column in matrix, which this column is CD4T+CD8T
sum_cd4_cd8 <- function(mat){
	new_mat <- matrix(NA, nrow(mat), ncol(mat)+1)
	new_mat[, 1:ncol(mat)] <- mat
	new_mat[, ncol(new_mat)] <- mat[, 3] + mat[, 4]	
	tmp_char <- union(colnames(mat),"cd4+8")
	colnames(new_mat) <- tmp_char
	return(new_mat)
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



