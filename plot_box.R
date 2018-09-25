library(ggplot2)

setwd("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/TIMER_rmse_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/EPIC_rmse_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/CIBER_rmse_list.RData")

load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/TIMER_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/EPIC_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/CIBER_333_list.RData")
load("C:/Users/wnchang/Documents/F/PhD_Research/2018_08_23_deconvolution_score/ICTD_333_list.RData")




cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
cancer_lib <- c("BRCA","COAD")


# 1 B
p_data <- c()
pdf("testB.pdf")
for(k in 1:length(cancer_lib)){
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

	itmp <- ICTD_rmse[[k]][["B_mark"]]
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


	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2","gold"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":B cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5))
    print(pp)

}
dev.off()



# 2 CD4T
p_data <- c()
pdf("testCD4T.pdf")
for(k in 1:length(cancer_lib)){
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

	itmp <- ICTD_rmse[[k]][["CD4_mark"]]
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


	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2","gold","grey"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":CD4T cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()

# 3 CD8T
p_data <- c()
pdf("testCD8T.pdf")
for(k in 1:length(cancer_lib)){
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

		itmp <- ICTD_rmse[[k]][["CD8_mark"]]
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


	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid","palegreen3","navajowhite2"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":CD8T cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()

# 4 Neutrophil
p_data <- c()
pdf("testNeutrophil.pdf")
for(k in 1:length(cancer_lib)){
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["Netro_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["Neutrophils"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)


	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		#geom_boxplot(fill= c("lightsteelblue","orchid"), alpha=0.8) + 
   		geom_violin(aes(fill = as.factor(algo)),show.legend=F) +
    	xlab("algo") + ggtitle(paste(cancer_str,":Neutrophil cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()


# 5 Macrophage
p_data <- c()
pdf("testMacrophage.pdf")
for(k in 1:length(cancer_lib)){
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

		itmp <- ICTD_rmse[[k]][["Macro_mark"]]
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


	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		#geom_boxplot(fill= c("lightsteelblue","orchid", "palegreen3", "navajowhite2", "gold","grey"), alpha=0.8) + 
   		geom_violin(aes(fill = as.factor(algo)),show.legend=F) +
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":Macrophage cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()


# 6 Dendritic cell
p_data <- c()
pdf("testDC.pdf")
for(k in 1:length(cancer_lib)){
	p_data <- c()
	cancer_str <- cancer_lib[k]
	ttmp <- TIMER_rmse[[k]][["DC_mark"]]
	rmse <- matrix(ttmp, length(ttmp), 1)
	method <- rep("timer", length(ttmp))
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(ttmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)


	ctmp <- CIBER_rmse[[k]][["Dendritic.cells.resting"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_DC_rest", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["Dendritic.cells.activated"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_DC_act", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)

	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid", "palegreen3"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) +
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":DC cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()

#7 NK 
p_data <- c()
pdf("testNK.pdf")
for(k in 1:length(cancer_lib)){
	p_data <- c()
	cancer_str <- cancer_lib[k]
	
	etmp <- EPIC_rmse[[k]][["NK_mark"]]
	rmse  <- matrix(etmp, length(etmp), 1)
	method <- rep("epic", length(etmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(etmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	

		itmp <- ICTD_rmse[[k]][["NK_mark"]]
	rmse  <- matrix(itmp, length(itmp), 1)
	method <- rep("ictd", length(itmp), 1)
	algo <- matrix(method, length(method), 1)
	rownames(rmse) <- names(itmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)	

	ctmp <- CIBER_rmse[[k]][["NK.cells.resting"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_NK_rest", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	ctmp <- CIBER_rmse[[k]][["NK.cells.activated"]]
	rmse <- matrix(ctmp, length(ctmp), 1)
	method <- rep("ciber_NK_act", length(ctmp), 1)
	algo <- matrix(method, length(ctmp), 1)
	rownames(rmse) <- names(ctmp)
	rmse <- cbind.data.frame(rmse, algo)
	p_data <- rbind(p_data, rmse)

	p_data <- as.data.frame(p_data)

	pp <- ggplot(p_data, aes(x=as.factor(algo), y=rmse)) + 
   		geom_boxplot(fill= c("lightsteelblue","orchid", "palegreen3", "navajowhite2"), alpha=0.8) + 
   		#geom_violin(aes(fill = as.factor(algo)),show.legend=F) + 
    	xlab("algo") + ylab("R2") + ggtitle(paste(cancer_str,":NK cell R2 comparison",sep="") ) +
    	#ylim(0, 0.2) +
    	theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    print(pp)

}
dev.off()

