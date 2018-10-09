rm(list = ls())
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

#choose just one of below dataset to run

# 1. GSE72056 melanoma
dataSet = "GSE72056"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE72056_tg_data_list.RData")
tg_data_list <- GSE72056_tg_data_list 
Cell_Prop <- Cell_Prop_GSE72056

# 4. GSE103322 HNC
dataSet <- "GSE103322"
list_flag <- T
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_06_28_singleCellSimulation/GSE103322_tg_data_list.RData")
tg_data_list <- GSE103322_tg_data_list 
Cell_Prop <- Cell_Prop_GSE103322





# set parameter
d=1		#we just choose one type of infiltratin data
if(list_flag == T)	bulk <- tg_data_list[[d]][[1]]
if(list_flag == F)	bulk <- tg_data_list[[d]]
if(list_flag == T) true_p <- tg_data_list[[d]][[2]]
if(list_flag == F) true_p <- Cell_Prop[[d]]	#??
colnames(bulk) <- paste("ss",c(1:ncol(bulk)), sep="" )
colnames(true_p) <- paste("ss", c(1:ncol(true_p)), sep="")

###################################
rm_zero_row <- function(bulk){

	row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == T)
	print(length(Zero))
	if(length(Zero) == 0) bulk <- bulk
	if(length(Zero) > 0)  bulk <- bulk[-Zero,]

	return(bulk)
}
############################
bulk <- rm_zero_row(bulk)

#prepare high correlation data
ccc <- cor(t(bulk), t(true_p))
dim(ccc)
colnames(ccc)

#######################################
find_correlated_marker <- function(ccc, cut){
	
	for(i in 1:ncol(ccc)){
		vv = ccc[which(ccc[,i] > cut ), i]
		vv_name = names(vv)

		assign(paste(colnames(ccc)[i],"_mark", sep=""), vv_name)
	}
	llll = list()
	for(i in 1:ncol(ccc)){
		llll[[i]] = get( paste(colnames(ccc)[i],  "_mark",  sep="") )
	}
	names(llll) = colnames(ccc)
	return(llll)

} 
##########################################

cutoff <- 0.9
cell_marker <- find_correlated_marker(ccc, cutoff)

#generate gene expression data XXX
XXX <- c()
for(i in 1:length(cell_marker)){
	XXX <- rbind(XXX, bulk[cell_marker[[i]], ])
}

#generate indicator matrix C
C <- matrix(0, nrow(XXX), length(cell_marker))
rownames(C) <- rownames(XXX)
colnames(C) <- names(cell_marker)
up <- 0
down <- 0
for(i in 1:length(cell_marker)){
	up <- down+1
	down <- up + length(cell_marker[[i]]) -1
	C[up:down, i] <- rep(1, length(cell_marker[[i]]))
}

bulk = XXX
ss_signature = C
save(bulk, ss_signature, true_p, file="103322_X_C.Rdata")

#####################################################################
#	Here, bulk and ss_signature are input for NMF.
#	However, we need to choose the cell type (column of bulk and 
#   column of ss_signature) we want to use.
#####################################################################


