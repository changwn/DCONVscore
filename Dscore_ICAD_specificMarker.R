rm(list = ls())

library(EPIC)
library("caret")
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

ICTD_rmse <- list()

#remove LAML
cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
cancer_lib <- c("BRCA","COAD","BRCA_TNBC")

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
	cell_marker <- extract_marker(signature,method="EPIC")
	ss_signature <- unique_signature(signature, 0.05)

	row_sub <- apply(ss_signature, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) ss_signature <- ss_signature
	if(length(Zero) > 0)  ss_signature <- ss_signature[-Zero,]
	ss_gene <- rownames(ss_signature)


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

	ICTD_rmse[[k]] <- cell_rmse
	names(ICTD_rmse)[k] <- cancer_str

}	# for

save(ICTD_rmse, file = "ICTD_333_list.RData")


