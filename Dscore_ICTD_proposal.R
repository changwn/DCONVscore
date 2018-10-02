
#COAD
#names(R1_filter_step1_results_new[[4]])[c(1,6,8,9,10)]

#TNBC
#names(R1_filter_step1_results_new[[4]])[c(2,3,4,6,9)]

#BRCA
names(R1_filter_step1_results_new[[4]])[c(1,3,5,7,8)]

ICTD_rmse <- list()
#cancer_str <- "BRCA"
#cancer_str <- "BRCA_TNBC"
cancer_str <- "COAD"
print(cancer_str)

#rm(list = ls())
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_10_01_data/TCGA-COAD_newpipeline201809.RData")

setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")



file_str <- paste("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA-", cancer_str, "_FPKM_T.RData", sep="")
load(file_str)

data_ttt <- filter_gene_name(data_t)
bulk <- data_ttt
# just use gene which can produce protein
bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
# remove whole zero row
bulk <- rm_zero_row(bulk)


# prepare marker
#dlll <- R1_filter_step1_results_new[[4]][c(1,3,5,7,8)]
dlll <- R1_filter_step1_results_new[[4]][c(1,6,8,9,10)]	#coad

if(F){
adipocytes_mark <- dlll[[1]]
B_mark <- dlll[[2]]
Myeloid_mark <- dlll[[3]]
Fibro_mark <- dlll[[4]]
TNK_mark <- dlll[[5]]

mark5 <- c(adipocytes_mark, B_mark, Myeloid_mark, Fibro_mark, TNK_mark)
XXX <- bulk[mark5, ]
C <- matrix(0, nrow(XXX), 5)
rownames(C) <- rownames(XXX)
colnames(C) <- c("adipocyte","B","Myeloid","Fibro","TNBCK")
length(adipocytes_mark)
length(B_mark)
length(Myeloid_mark)
length(Fibro_mark)
length(TNK_mark)
C[1:17, 1] <- rep(1, 17)
C[18:35, 2] <- rep(1, 18)
C[36:51, 3] <- rep(1, 16)
C[52:70, 4] <- rep(1, 19)
C[71:90, 5] <- rep(1, 20)
}
if(T){
Fibro_mark <- dlll[[1]]
B_mark <- dlll[[2]]
T_mark <- dlll[[3]]
Myeloid_mark <- dlll[[4]]
Endothelial_mark <- dlll[[5]]

mark5 <- c(Fibro_mark, B_mark, T_mark, Myeloid_mark, Endothelial_mark)
XXX <- bulk[mark5, ]
C <- matrix(0, nrow(XXX), 5)
rownames(C) <- rownames(XXX)
colnames(C) <- c("Fibro","B","T","Myeloid","Endo")
length(Fibro_mark)
length(B_mark)
length(T_mark)
length(Myeloid_mark)
length(Endothelial_mark)
C[1:7, 1] <- rep(1, 7)
C[8:24, 2] <- rep(1, 17)
C[25:44, 3] <- rep(1, 20)
C[45:55, 4] <- rep(1, 11)
C[56:65, 5] <- rep(1, 10)	
}

bulk = XXX
ss_signature = C

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

save(ICTD_rmse, file = "ICTD_333_list.RData")






















########################## function


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


rm_zero_row <- function(bulk){

	row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) bulk <- bulk
	if(length(Zero) > 0)  bulk <- bulk[-Zero,]

	return(bulk)
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


 