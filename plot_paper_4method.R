library(ggplot2)
library(gridExtra)
#library(plotly)


setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")


load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/TIMER_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/EPIC_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/CIBER_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/ICTD_333_list.RData")


load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/TIMER_103322_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/EPIC_103322_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/CIBER_103322_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/ICTD_103322_list.RData")



cancer_lib <- c("COAD","BRCA")
#cancer_lib <- c("GSE103322")
#cancer_str <- c("103322")
pp <- list()
count <- 1

k=1


# 1 B
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["B_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	etmp <- EPIC_rmse[[k]][["B_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	itmp <- ICTD_rmse[[k]][["B"]]
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["B.cells.naive"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_B_naive", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["B.cells.memory"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_B_memory", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)


	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2","gold"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":B cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),,axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())

    count <- count + 1

# 2 CD4T
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["CD4_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	etmp <- EPIC_rmse[[k]][["CD4_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	itmp <- ICTD_rmse[[k]][["TNBCK"]]
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)


	ctmp <- CIBER_rmse[[k]][["T.cells.CD4.naive"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_cd4_naive", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["T.cells.CD4.memory.resting"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_cd4_mem_rest", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["T.cells.CD4.memory.activated"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_cd4_mem_act", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)


	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2","gold","grey"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":CD4T cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())
    count = count + 1

# 3 CD8T
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["CD8_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	etmp <- EPIC_rmse[[k]][["CD8_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

		itmp <- ICTD_rmse[[k]][["TNBCK"]]
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["T.cells.CD8"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_cd8", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)

	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":CD8T cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())
    	count = count + 1

# 4 Macrophage
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["Macro_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	etmp <- EPIC_rmse[[k]][["Macro_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

		itmp <- ICTD_rmse[[k]][["Myeloid"]]	###myeloid
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	

	ctmp <- CIBER_rmse[[k]][["Macrophages.M0"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_M0", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["Macrophages.M1"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_M1", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["Macrophages.M2"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_M2", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)


	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid", "palegreen3", "navajowhite2", "gold","grey"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) +
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":Macrophage cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())
    	count = count + 1

#5 Fibroblat 
	p_data <- c()
	cancer_str <- cancer_lib[k]
	
	etmp <- EPIC_rmse[[k]][["CAFs_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	

		itmp <- ICTD_rmse[[k]][["Fibro"]]
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	


	p_data <- as.data.frame(p_data)

	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("orchid", "palegreen3"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":Fibroblast cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())
    count = count + 1

#6 Endothelial
	p_data <- c()
	cancer_str <- cancer_lib[k]
	
	etmp <- EPIC_rmse[[k]][["endoth_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	

		itmp <- ICTD_rmse[[k]][["adipocyte"]]	###adipocyte
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	


	p_data <- as.data.frame(p_data)

	pp[[count]] <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("orchid", "palegreen3"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":Endothelial cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(), 
    												axis.text.x=element_blank(),
                      								axis.ticks.x=element_blank())
    count = count + 1



#pdf("1111.pdf", width = 30, height = 12)
pdf("R2_COAD_4method.pdf", width = 24, height = 5)
do.call(grid.arrange, c(pp,nrow=1 ) )
dev.off()
















