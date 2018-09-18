#using simulate single cell data to verify TIMER pipeline is ture
#CANNOT verify because it just work for TCGA data


load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_28_Brianna/signature_matrix_from_three_tools.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA_ensem_annotation.RData")

##############################################
file_str <- "C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA-COAD_FPKM_T.RData"
load(file_str)

data_ttt <- filter_gene_name(data_t)
bulk <- data_ttt

# just use gene which can produce protein
bulk <- bulk[intersect(rownames(bulk),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]

# remove whole zero row
row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk <- bulk
if(length(Zero) > 0)  bulk <- bulk[-Zero,]

#com <- intersect(rownames(bulk), gene_name)
#bulk <- bulk[com,]
bulk <- rbind.data.frame(colnames(bulk), bulk)
bulk <- cbind.data.frame( rownames(bulk) , bulk  )
rownames(bulk) <- c()
colnames(bulk) <- c()

write.csv(bulk, file="coad.csv", row.names = F, col.names = F)

serve ERROR



##########################################
bulk <- read.csv("./exampleForLUAD.csv", header = T)
rownames(bulk) <- bulk[, 1]
bulk <- bulk[, -1]

cancer_str = "LUAD"
cancer_str <- tolower(cancer_str)
print(cancer_str)

dd = bulk

Error in cor.test.default(tmp.dd[x, ], as.numeric(AGP[colnames(tmp.dd),  : 
  not enough finite observations


#######################################################
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
Prop_t <- Cell_Prop[[d]]

# delete 0 row
row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
Zero <- which(row_sub == T)
if(length(Zero) == 0) bulk1 <- bulk
if(length(Zero) > 0)  bulk1 <- bulk[-Zero,]

cancer_str = "SKCM"
cancer_str <- tolower(cancer_str)
print(cancer_str)

timer_result <- TIMER(cancer_str)

Error in cor.test.default(tmp.dd[x, ], as.numeric(AGP[colnames(tmp.dd),  : 
  not enough finite observations


#names(timer_result)
#dim(timer_result[[1]])
Prop_TIMER <- timer_result[[1]]

cor(Prop_TIMER, t(Prop_t) )
