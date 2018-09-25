library(ggplot2)
library(gridExtra)

setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/TIMER_rmse_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/EPIC_rmse_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/CIBER_rmse_list.RData")


p_data <- c()
pp <- list()

cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
#cancer_lib <- c("BLCA","COAD")

#B
for(k in 1:length(cancer_lib)){
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


	pp[[k]]<- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_violin(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2"), alpha=0.8) + 
    	xlab("algo") + ggtitle(paste(cancer_str,":B cell RMSE comparison",sep="") ) +
    	ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5))
    
}

#pdf("1111.pdf", width = 30, height = 12)
do.call(grid.arrange, c(pp,nrow=4 ) )
#dev.off()
