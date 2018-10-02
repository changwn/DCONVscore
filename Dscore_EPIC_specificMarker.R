rm(list = ls())

library(EPIC)
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

EPIC_rmse <- list()

#remove LAML
cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
cancer_lib <- c("COAD","BRCA_TNBC")

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
	bulk <- rm_zero_row(bulk)


#----------------------
# seed <- MR? 98 genes
# pool <- high correlation
#----------------------

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
	ss_signature <- unique_signature(signature)
	ss_signature <- rm_zero_row(ss_signature)
	bulk1 <- bulk1[rownames(ss_signature), ]

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
	#bulk_est222 = bulk_est
	#}
	#3{
	#CC = do_correlation(Prop_EPIC[,1:7], bulk1, ss_signature)
	SS <- do_regression_unique(Prop_EPIC, bulk1, ss_signature)
	bulk_est <- SS %*% t(Prop_EPIC)
	#bulk_est333 = bulk_est
	#}



	#a. direct rmse
	#rmse_gene <- RMSE_two_mat(bulk1, bulk_est)
	#b.  rmse of radio
	#{
	#rmse_gene <- RMSE_two_mat(bulk1, bulk_est)
	#gene_vira <- row_vira(bulk1)
	#rmse_gene <- rmse_gene / gene_vira
	#}
	# c. R2
	rmse_gene <- R2_two_mat(bulk_est, bulk1)

	


	####debug
	#elist <- list(bulk1, signature[rownames(ss_signature),], Prop_EPIC, ss_signature)
	
	#extract marker corresponding rmse
	cell_rmse <- list()
	for(i in 1:length(cell_marker)){
		cell_rmse[[i]] <- rmse_gene[cell_marker[[i]], ]
	}
	names(cell_rmse) <- names(cell_marker)
	#plot
	#ggred <- 
	#ggblue <- 
	#pdf_str <- paste(cancer_str, "_epic_cell_rmse_unique.pdf", sep="")
	#pdf_str <- "Timer_cell_rmse.pdf"
	#pdf(file = pdf_str)
	flag_plot <- F
	if(flag_plot == T){
	pdf_str <- paste(cancer_str, "_epic_cell_rmse_unique.pdf", sep="")
	pdf(file = pdf_str)	
	for(i in 1:length(cell_marker)){
		str <- paste(cancer_str, names(cell_rmse)[i], "rmse_unique", sep="-")
		pp <- barplot(cell_rmse[[i]], main=str, col = c("lightblue"), names.arg="")
		text(pp, -0.002, srt = 45, adj= 1, xpd = TRUE, labels = names(cell_rmse[[i]]) , cex=0.5)
	}
	dev.off()
	}


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


	EPIC_rmse[[k]] <- cell_rmse
	names(EPIC_rmse)[k] <- cancer_str

}	# for

save(EPIC_rmse, file = "EPIC_333_list.RData")






#gene_Smean <- mean(rmse_gene)		#init mean score
#gene_Smedian <- median(rmse_gene)	#init median score











#-----------------------------
#function

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

do_regression <- function(proportion=Prop_EPIC, data=bulk1){
	P_2nd <- c()
	n_gene <- nrow(data)
	g_sample <- ncol(data)
	for(i in 1:n_gene){
		coeff <- getFractions.Abbas(proportion, t(data)[, i])
		P_2nd <- rbind(P_2nd, coeff)
	}

	return(P_2nd)
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
		if(length(cell) > 1) stop("return 2 cell type!")
	 	pp <- proportion[, cell]
		y <- data[i, ]
		ddd <- cbind(pp, y)
		ddd <- as.data.frame(ddd)
		coeff <- lm(y~pp + 0, ddd)	#NOTE : the order is matter!
		P_2nd[i, cell] <- coeff[[1]]
	}	
	return(P_2nd)
 }

do_correlation <- function(proportion, data, indicate){
	n_gene <- nrow(data)
	n_sample <- ncol(data)	
	for(i in 1:n_gene){
		cell <- which(indicate[i, ] == 1)
		if(length(cell) > 1) stop("return 2 cell type!")
	 	pp <- proportion[, cell]
	 	y <- data[i, ]
	 	print(cor(pp,y) )


	}

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

row_vira <- function(aaa){
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		ave <- mean(aaa[i, ])
		vc[i] <- sqrt(sum( (aaa[i, ] - ave)^2 ) )

	}
	return(vc)
}


extract_marker <- function(aaa, method="EPIC"){
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


