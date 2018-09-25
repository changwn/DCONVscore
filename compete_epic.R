rm(list = ls())

library(EPIC)
library("caret")
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

EPIC_rmse <- list()

#remove LAML
cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
cancer_lib <- c("BLCA","COAD")

for(k in 1:length(cancer_lib)){
	cancer_str <- cancer_lib[k]
	print(cancer_str)

	file_str <- paste("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA-", cancer_str, "_FPKM_T.RData", sep="")
	load(file_str)

	data_ttt <- filter_gene_name(data_t)
	bulk <- data_ttt

	# just use gene which can produce protein
	bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

	# remove whole zero row
	row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) bulk <- bulk
	if(length(Zero) > 0)  bulk <- bulk[-Zero,]

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
	ss_signature <- unique_signature(signature, 0.05)

	row_sub <- apply(ss_signature, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) ss_signature <- ss_signature
	if(length(Zero) > 0)  ss_signature <- ss_signature[-Zero,]
	ss_gene <- rownames(ss_signature)

	bulk1 <- bulk[ss_gene,]

	#colnames(bulk1) <- paste("ssample",c(1:ncol(bulk1)), sep="")
	out1 <- EPIC(bulk1)
	print("finish EPIC")
	Prop_EPIC <- out1$cellFraction[,1:7]

	# (1) use whole prop or (2) use cell specific prop	
	#1{
	#S <- do_regression(Prop_EPIC, bulk1)		#!: use small data do regression
	#bulk_est <- S %*% t(Prop_EPIC)
	#}
	#2{
	#bulk_est <- ss_signature %*% t(Prop_EPIC[,1:7]) 
	#rownames(bulk_est) <- rownames(bulk1)
	#}
	#3{
	Sc <- do_regression_unique(Prop_EPIC, bulk1, ss_signature)
	bulk_est <- Sc %*% t(Prop_EPIC)	
	#}	

	#1 direct rmse
	#rmse_gene <- RMSE_two_mat(bulk1, bulk_est)
	#2 rmse of radio
	#{
	#rmse_gene <- RMSE_two_mat(bulk1, bulk_est)
	#gene_vira <- row_vira(bulk1)
	#rmse_gene <- rmse_gene / gene_vira
	#}
	rmse_gene <- R2_two_mat(bulk_est, bulk1)

	#extract marker corresponding rmse
	cell_rmse <- list()
	for(i in 1:length(cell_marker)){
		cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
	}
	names(cell_rmse) <- names(cell_marker)

	flag_plot <- T
if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_epic_cell_R2_unique.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "R2", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()

}	


############### NMF
	data_t = log2(bulk + 1)
	NMF_indi_all = ss_signature

X1=data_t[match(rownames(NMF_indi_all),rownames(data_t)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
###########

###########Parameter settings
theta=0.5 ##penalty parameter for constraints on NMF_indi_all
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

	rmse_gene <- R2_two_mat(U %*% t(V), X1)

	#extract marker corresponding rmse
	cell_rmse <- list()
	for(i in 1:length(cell_marker)){
		cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
	}
	names(cell_rmse) <- names(cell_marker)

	flag_plot <- T
	if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_ICAD_cell_R2_unique.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "R2", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()
	}



}	# for




















#--------------------------------------------------------


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


unique_signature <- function(aaa, cutoff = 0.05){

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


RMSE_one_vector <- function(a){
	rmse <- sqrt(sum(a^2)/length(a))
	return(rmse)
}

RMSE_two_vector <- function(a, b){
	diff <- a - b
	rmse <- sqrt(sum(diff^2)/length(diff))
	return(rmse)
}

library("caret")
R2_two_vector <- function(predict, actual){
	R2 <- caret::postResample(predict, actual)[["Rsquared"]]
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

