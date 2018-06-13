getEstimation<-function(N){
  
  directory <-paste0("D:\\insa\\inserm\\Simulation\\result\\md_",1,"\\projectSimulSimon")
  setwd(directory)
  Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
  datas<-Data[,2]
  
  for(i in 2:N){
    tryCatch({
      directory <-paste0("D:\\insa\\inserm\\Simulation\\result\\md_",i,"\\projectSimulSimon")
      setwd(directory)
      Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
      datas<-cbind(datas,Data[,2])
    },error=function(e){
      print(i)
    }
    )
    
  }
  return(datas)
}

gg<-getEstimation(500)
dim(gg)
vecName<-c("¦Ì_VA","¦Ì_k","¦Ì_be","¦Ì_Emax","¦Ì_ED50","w2_VA","w2_k","w2_be","w2_Emax","w2_ED50","sigma_i")

SEemp<-c()
for(i in 1:11){
  SEemp <- c(SEemp,sqrt(var(gg[i,])))
  # boxplot(gg[i,])
}

names(SEemp) <- vecName
# SEemp[1:5]/c(55,0.005,0.2,30,150)*100

paramValue<-c(55,0.005,0.2,30,150)
omega<-c(0.07,0.5,1,150,0)
sigma_i<-sqrt(28)
ree<-(gg[1,]-paramValue[1])/paramValue[1]*100
for(i in 2:5){
  ree<-cbind(ree,(gg[i,]-paramValue[i])/paramValue[i]*100)
}
for(i in 1:5){
  ree<-cbind(ree,(gg[5+i,]-omega[i])/omega[i]*100)
}
ree<-cbind(ree,(gg[11,]-sigma_i)/sigma_i*100)

names(ree)<-vecName
#Boxplots ree
library(Rlab)
bplot(ree[1,],style="quantile",labels=c("¦Ì_VA","¦Ì_k","¦Ì_be","¦Ì_Emax","¦Ì_ED50","w2_VA","w2_k","w2_be","w2_Emax","w2_ED50","sigma_i"))
abline(h=0,col='red')

RB<-apply(ree,2,mean)
names(RB)<-vecName

RSE<-100 * SEemp / c(paramValue,omega,sigma_i)[which(c(paramValue,omega,sigma_i)!=0)]
RSE<-RSE[-c(10)]


source("D:\\insa\\inserm\\FIMmixte.r")
fimVA<-funFIMem("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,168,644,672),c(0,7,644,672),c(0,168,196,672)),c(2,2,2,1,2),c(0,100,1000),c(0.2,0.5,0.3),300)
RSEfim<-fimVA[[5]]

RSE2<-rbind(RSE,RSEfim)
barplot(RSE2,beside = T,legend=T)

par(mfrow=c(2,2))
for(i in 1:11){
  hist(gg[i,],main=paste0("histogram de ",vecName[i]))
}

# 
# boxplot(ylab=list("REE(%)",cex=1.4),main=list("VA",cex=2),ree,las=1, 
#         names=c("¦Ì_VA","¦Ì_k","¦Ì_be","¦Ì_Emax","¦Ì_ED50","w2_VA","w2_k","w2_be","w2_Emax","w2_ED50","sigma_i"))
