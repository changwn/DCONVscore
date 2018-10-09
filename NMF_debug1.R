ss_signature
bulk
rownames(bulk)
rownames(ss_signature)
hist(log2(bulk+1))

hist(data_t)

library(gplots)

colors = c(-100:100)/100
my_palette <- colorRampPalette(c("red","white", "blue"))(n =200)

aaa<-cor(t(data_t))
heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)


kk<-5
tg_ids<-c()
for(i in 1:ncol(ss_signature))
{
	tg_ids0<-which(ss_signature[,i]==1)
	tg_ids0<-tg_ids0[1:min(kk,length(tg_ids0))]
	tg_ids<-c(tg_ids,tg_ids0)
}


ss1<-ss_signature[tg_ids,]
dd1<-data_t[tg_ids,]

aaa<-cor(t(dd1))
heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)

############### NMF
NMF_indi_all = ss1	#NMF input 2

X1=dd1[match(rownames(NMF_indi_all),rownames(dd1)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
###########

###########Parameter settings
theta=500 ##penalty parameter for constraints on NMF_indi_all
indiS_method="nonprdescent" ##the updating scheme for the structural constraints
iter=2000
alpha=beta=gamma=roh=0
roh=gamma=0.0001
UM=VM=NULL
qq=1
epslog=6
nPerm=4
initial_U=initial_V=NULL
mscale=1
###########

ttt1=qnmf_indisS_all_revise(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale)

library(NMF)
res<-nmf(X1,rank=5,method="snmf/l",nrun=10)
basismap(res,Rowv=NA,Colv=NA)
coefmap(res,Rowv=F)
 basis(res)
coef(res)

plot(ttt1$U[,5])

plot(ttt1$U[,1])
ttt1$V

