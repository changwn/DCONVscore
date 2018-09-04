
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
# 2. GSE70630 oligodendroglioma % in brain?
dataSet = "GSE70630"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE70630_tg_data_list.RData")
tg_data_list <- GSE70630_tg_data_list 
Cell_Prop <- Cell_Prop_GSE70630
# 3. GSE89567 astrocytoma
dataSet <- "GSE89567"
list_flag <- F
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE89567_tg_data_list.RData")
tg_data_list <- GSE89567_tg_data_list 
Cell_Prop <- Cell_Prop_GSE89567
# 4. GSE103322 HNC
dataSet <- "GSE103322"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE103322_tg_data_list.RData")
tg_data_list <- GSE103322_tg_data_list 
Cell_Prop <- Cell_Prop_GSE103322


# 5. GSE75688 breast cancer
dataSet <- "GSE75688"
list_flag <- F
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_29_singleCellSimulated/GSE75688_tg_data_list.RData")
tg_data_list <- GSE75688_tg_data_list 
Cell_Prop <- Cell_Prop_GSE75688
# 6. GSE81861 colorectal tumors
dataSet <- "GSE81861"
list_flag <- F
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_29_singleCellSimulated/GSE81861_tg_data_list.RData")
tg_data_list <- GSE81861_tg_data_list 
Cell_Prop <- Cell_Prop_GSE81861


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
print(Dscore_whole)

# get TOP 10: sort the gene based on correlation and then extract top 10 gene
dd_copy_order <- sort(dd_copy, decreasing = T)
base <- names(dd_copy_order[1:10])
out_base <- EPIC(bulk1, sigGenes = base)
prop_base <- out_base$cellFraction[, 1:7]
Dscore_init <- cal_Dscore_abs(prop_base, bulk1)
print(Dscore_init)

# Iteration: test remain gene 
remainSet <- setdiff(names(dd_copy_order), base)
top_gene_add <- base
if(dataSet == "GSE72056") Dscore_mat <- matrix(NA, length(remainSet), 8+2)
if(dataSet == "GSE70630") Dscore_mat <- matrix(NA, length(remainSet), 1)
if(dataSet == "GSE89567") Dscore_mat <- matrix(NA, length(remainSet), 1)
if(dataSet == "GSE103322") Dscore_mat <- matrix(NA, length(remainSet), 9+2)
if(dataSet == "GSE75688") Dscore_mat <- matrix(NA, length(remainSet), 7+2)
if(dataSet == "GSE81861") Dscore_mat <- matrix(NA, length(remainSet), 7+2)
#Dscore_mat <- c()
increase_gene <- c()
for(i in 1:length(remainSet)){
	print(i)
	top_gene_add <- union(top_gene_add, remainSet[i])
	out_ep <- EPIC(bulk1, sigGenes=top_gene_add)
	prop_add <- out_ep$cellFraction[,1:8]
	prop_add <- sum_cd4_cd8(prop_add)
	prop_true <- t(Cell_Prop[[d]])
	corr_add <- cor(prop_add, prop_true)
	if(dataSet == "GSE72056"){	
		Dscore_mat[i, 1] <- corr_add[1, 1]	#B cell
		Dscore_mat[i, 2] <- corr_add[9, 2]	#T
		Dscore_mat[i, 3] <- corr_add[8, 3]	#tumor
		Dscore_mat[i, 4] <- corr_add[6, 4]	#macrophage
		Dscore_mat[i, 5] <- corr_add[2, 5]	#Fibroblast
		Dscore_mat[i, 6] <- corr_add[5, 6]	#endothelial
		Dscore_mat[i, 7] <- mean(Dscore_mat[i, 1:6])	#average
		tmp <- (sum(Dscore_mat[i,1:6]) - Dscore_mat[i, 3]) / 5
		Dscore_mat[i, 8] <- tmp  			#ave-5
		Dscore_mat[i, 9] <- corr_add[3, 2]	#cd4t
		Dscore_mat[i,10] <- corr_add[4, 2]	#cd8t
		if(i==1) colnames(Dscore_mat) <- c("B","T","tumor","macrophage","Fibroblast","endothelial","ave_all","ave_5","cd4t","cd8t")

	}
	if(dataSet == "GSE70630"){
		Dscore_mat[i, 1] <- corr_add[6, 2]
		colnames(Dscore_mat) <- c("macrophages_macroglio")
	}
	if(dataSet == "GSE89567"){
		Dscore_mat[i, 1] <- corr_add[6, 5]
		colnames(Dscore_mat) <- c("macroglio")	
	}
	if(dataSet == "GSE103322"){
		Dscore_mat[i, 1] <- corr_add[8 ,1]	#tumor
		Dscore_mat[i, 2] <- corr_add[9, 2]	#CD4+8
	
		Dscore_mat[i, 3] <- corr_add[1, 3]	#B
		Dscore_mat[i, 4] <- corr_add[2, 5]	#fibroblast
		Dscore_mat[i, 5] <- corr_add[5, 6]	#endothelial
		Dscore_mat[i, 6] <- corr_add[5, 7]	#epithelial?
		Dscore_mat[i, 7] <- corr_add[6, 8]	#macrophage
		Dscore_mat[i, 8] <- mean(Dscore_mat[i, 1:7])	#ave_all
		tmp <- (sum(Dscore_mat[i,1:7]) - Dscore_mat[i, 1]) / 6	#ave_6
		Dscore_mat[i,9] <- tmp
		Dscore_mat[i,10]<- corr_add[3, 2]	#cd4t
		Dscore_mat[i,11]<- corr_add[4, 2]	#cd8t
		if(i==1) colnames(Dscore_mat) <- c("tumor","cd4+8t","B","fibroblast","endothelial","epithelial","macrophage","ave_all","ave_6","cd4t","cd8t")
	}
	if(dataSet == "GSE75688"){
		Dscore_mat[i, 1] <- corr_add[1, 1]	#B
		Dscore_mat[i, 2] <- corr_add[9, 2]	#CD4+8
		
		Dscore_mat[i, 3] <- corr_add[8, 3]	#tumor
		Dscore_mat[i, 4] <- corr_add[6, 4]	#myeloid_macrophage
		Dscore_mat[i, 5] <- corr_add[2, 5]	#stromal(fibroblast)
		Dscore_mat[i, 6] <- (sum(Dscore_mat[i,1:5]) - Dscore_mat[i, 3]) / 4
		Dscore_mat[i, 7] <- mean(Dscore_mat[i, 1:5])
		Dscore_mat[i, 8] <- corr_add[3, 2]	#cd4t
		Dscore_mat[i, 9] <- corr_add[4, 2]	#co8t
		if(i==1) colnames(Dscore_mat) <- c("B","cd4+8t","tumor","myeloid_macrophage","fibroblast","ave_5","ave_all","cd4t","cd8t")
	}
	if(dataSet == "GSE81861"){
		Dscore_mat[i, 1] <- corr_add[9, 1]	#cd4+8
		
		Dscore_mat[i, 2] <- corr_add[5, 2]	#epithelial_endothelial
		Dscore_mat[i, 3] <- corr_add[1, 3]	#B
		Dscore_mat[i, 4] <- corr_add[6, 5]	#macrophage
		Dscore_mat[i, 5] <- corr_add[2, 6]	#fibroblast
		Dscore_mat[i, 6] <- corr_add[5, 7]	#endothelial
		Dscore_mat[i, 7] <- mean(Dscore_mat[i, 1:6])
		Dscore_mat[i, 8] <- corr_add[3, 1]	#cd4t
		Dscore_mat[i, 9] <- corr_add[4, 1]	#cd8t
		if(i==1) colnames(Dscore_mat) <- c("cd4+8t","err:epithelial","B","macrophage","fibroblast","endothelial","ave_6","cd4t","cd8t")
	}
}

#plot
#par(mfrow = c(2,1))
#matplot(Dscore_mat, type=c("l"), pch = 1, col = 1:8)
#legend("bottomright", legend = colnames(Dscore_mat), col=1:8, pch=1)
if(dataSet == "GSE72056"){
	matplot(Dscore_mat[,-c(3,7)], type=c("l"), pch = 1, col = 1:8, main="GSE72056 melanoma, cor score")
	leg <- union(colnames(Dscore_mat)[1:2],colnames(Dscore_mat)[4:6])
	leg <- union(leg, colnames(Dscore_mat)[8:10])
	legend("bottomright", legend = leg, col=1:8, pch=1)
}
if(dataSet == "GSE70630"){
	plot(Dscore_mat, type = 'l', main = "GSE70630 oligodendroglioma")
	legend("bottomright", legend = colnames(Dscore_mat))
}
if(dataSet == "GSE89567"){
	plot(Dscore_mat, type = 'l', main = "GSE89567 astrocytoma")
	legend("bottomright", legend = colnames(Dscore_mat))
}
if(dataSet == "GSE103322"){
	#matplot(Dscore_mat, type='l', pch=1, col=1:10, main="GSE103322 HNC,cor score")
	#legend("bottomright", legend = colnames(Dscore_mat), col=1:10, pch=1)
	tmp_mat <- cbind(Dscore_mat[, 2:7], Dscore_mat[,9:11])
	matplot(tmp_mat, type='l', pch=1, col=1:9, main="GSE103322 HNC,cor score")
	leg <- union(colnames(Dscore_mat)[2:7],colnames(Dscore_mat)[9:11])
	legend("bottomright", legend = leg, col=1:9, pch=1)
}
if(dataSet == "GSE75688"){
	tmp_mat <- cbind(Dscore_mat[,1:2], Dscore_mat[,4:6])
	tmp_mat <- cbind(tmp_mat, Dscore_mat[,8:9])
	matplot(tmp_mat, type='l', pch=1, col=1:7, main="GSE75688 breast,cor score")
	leg <- union(colnames(Dscore_mat)[1:2], colnames(Dscore_mat)[4:6])
	leg <- union(leg, colnames(Dscore_mat)[8:9])
	legend("bottomright", legend = leg, col=1:7, pch=1)
}
if(dataSet == "GSE81861"){
	matplot(Dscore_mat, type="l", pch=1, col=1:7, main="GSE81861 colorectal,cor score")
	legend("bottomright", legend = colnames(Dscore_mat), col = 1:7, pch = 1)
}



#-----------------------------------------

# First, use explanation idea to give the order of all gene,
# then, use this order to calculate the correlation scoreby adding one gene each time.
# 1. if score not decend, means that above order not work?
# 2. the group of gene (row space) cannot be evaluate by adding one each time?
# 3. the cibersort list not work?

# question: 
# 1.how to give a initial set that can fully cover useful immune and cancer cell information,
# 2.hwo to give a sort method that let correlation score descend as gene number accumulated.

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



