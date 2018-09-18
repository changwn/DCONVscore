
library(EPIC)
setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")


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

bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

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

new_data <- bulk[commonSiga,]

# delete 0 row
row_sub <- apply(new_data, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk1 <- new_data
if(length(Zero) > 0)  bulk1 <- new_data[-Zero,]


colnames(bulk1) <- paste("ssample",c(1:ncol(bulk1)), sep="")
out1 <- EPIC(bulk1)
Prop_EPIC <- out1$cellFraction


ccc <- cor(Prop_EPIC, t(bulk1))
rownames(ccc)

#corr_score <- matrix(NA, nrow(ccc), ncol(ccc))
#for(i in 1:nrow(ccc)){
#	corr_score[i, ] <- sort(ccc[i, ], decreasing = T)
#}
#rownames(corr_score) <- rownames(ccc)

B_score <- sort(ccc[1, ], decreasing = T) 
F_score <- sort(ccc[2, ], decreasing = T) 
CD4_score <- sort(ccc[3, ], decreasing = T) 
CD8_score <- sort(ccc[4, ], decreasing = T) 
endo_score <- sort(ccc[5, ], decreasing = T) 
macroph_score <- sort(ccc[6, ], decreasing = T) 
NK_score <- sort(ccc[7, ], decreasing = T) 
tumor_score <- sort(ccc[8, ], decreasing = T) 

#intersect(names(B_score[1:10]), names(F_score[1:10]))
#intersect(names(B_score[1:20]), names(CD4_score[1:20]))
#Reduce(intersect, list(names(B_score[1:20]),  names(CD4_score[1:20]) ))

sets <- list(B_score, F_score, CD4_score, CD8_score, endo_score, macroph_score, NK_score)
Reduce(union, sets)

# correalation with  matrix V (in svd) of original data
svd_result <- svd(bulk1)
u <- svd_result$u
v <- svd_result$v
plot(svd(bulk1)$d)
dim(cor(Prop_EPIC, t(bulk1))  )
dim(cor(Prop_EPIC, v) )

# we assume the number of cell type will not exceed to 10, thus we just
# use top 10 basis. (e.g. first 10 column in V)
cor_with_basis <- cor(Prop_EPIC, v[, 1:10])
apply(cor_with_basis, 1, max)



############
#  test pca
###########


pca_sample <- prcomp(bulk1)
names(pca_sample)
dim(pca_sample$rotation)
dim(pca_sample$x)
summary(pca_sample)

pca_gene <- prcomp(t(bulk1))
dim(pca_gene$rotation)
dim(pca_gene$x)
summary(pca_gene)





