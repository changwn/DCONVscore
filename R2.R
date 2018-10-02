ls()

X<-clist[[1]]
S<-clist[[2]]
P<-clist[[3]]
C<-clist[[4]]

i<-1
cor(X[i,],S%*%P[i,])

dim(S)
dim(P)
#hist((S%*%t(P))[i,])

cor((S%*%t(P))[i,],X[i,])
ddd<-c()
for(i in 1:nrow(X))
{
	aaa<-summary(lm(X[i,]~(S%*%t(P))[i,]+0))
	ddd<-rbind(ddd,c(aaa[[4]][1,1],aaa[[8]]))
}

ddd2<-c()
for(i in 1:nrow(X))
{
	#tg_ids<-as.character(which(C[i,]==1))
tg_ids<-which(C[i,]==1)
	aaa<-summary( lm( X[i,]~t(  S[i,tg_ids] %*%t(P[,tg_ids])  )+0 )  )
	b <- R2_two_vector(S[i,tg_ids] %*%t(P[,tg_ids]),  X[i,])
	ddd2<-rbind(ddd2,c(aaa[[4]][1,1],aaa[[8]],b))
}



par(mfcol=c(1,3))
fff<-list()
ggg<-ddd
for(i in 1:ncol(C))
{
	fff[[i]]<-ggg[which(C[,i]==1),2]
}
names(fff)<-colnames(C)
boxplot(fff,las=2)

fff<-list()
ggg<-ddd2
for(i in 1:ncol(C))
{
	fff[[i]]<-ggg[which(C[,i]==1),2]
}
names(fff)<-colnames(C)
boxplot(fff,las=2)

####################################


X<-elist[[1]]
S<-elist[[2]]
P<-elist[[3]]
C<-elist[[4]]

dim(S)
dim(P)
#hist((S%*%t(P))[i,])

ddd2<-c()
for(i in 1:nrow(X))
{
	tg_ids<-colnames(C)[which(C[i,]==1)]
	aaa<-summary(lm(X[i,]~as.matrix(P[,tg_ids])+0))
	ddd2<-rbind(ddd2,c(aaa[[4]][1,1],aaa[[8]]))
}

fff<-list()
ggg<-ddd2
for(i in 1:ncol(C))
{
	fff[[i]]<-ggg[which(C[,i]==1),2]
}
names(fff)<-colnames(C)
boxplot(fff,las=2)



