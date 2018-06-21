getEstimation<-function(N){
  
  directory <-paste0("D:\\yuxin\\Simulation2\\resultCDC3\\md_",1,"\\project3")
  #directory <-paste0("C:\\yuxinCDC\\resultCDC2\\md_",54,"\\projectSimulSimon")
  setwd(directory)
  Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
  datas<-Data[,2]
  
  for(i in 2:N){
    tryCatch({
      directory <-paste0("D:\\yuxin\\Simulation2\\resultCDC3\\md_",i,"\\project3")
      #directory <-paste0("C:\\yuxinCDC\\resultCDC2\\md_",i,"\\projectSimulSimon")
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
gg<-gg[-c(10),]
vecName<-c("µ_VA0","µ_k","µ_be","µ_Emax","µ_ED50","w2_VA","w2_k","w2_be","w2_Emax","sigma_i")
#vecName<-c(expression(mu["VA"]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]^{2}),expression(omega["k"]^{2}),expression(omega[beta]^{2}),expression(omega["E"["max"]]^{2}),expression(sigma["i"]))
SEemp<-c()
for(i in 1:10){
  SEemp <- c(SEemp,sqrt(var(gg[i,])))
  # boxplot(gg[i,])
}

names(SEemp) <- vecName

paramValue<-c(55,0.005,0.2,30,150)
omega<-c(0.07,0.5,1,150)
sigma_i<-sqrt(28)
ree<-(gg[1,]-paramValue[1])/paramValue[1]*100
for(i in 2:5){
  ree<-cbind(ree,(gg[i,]-paramValue[i])/paramValue[i]*100)
}
for(i in 1:4){
  ree<-cbind(ree,(gg[5+i,]-omega[i])/omega[i]*100)
}
ree<-cbind(ree,(gg[10,]-sigma_i)/sigma_i*100)

names(ree)<-vecName
#Boxplots ree
#library(Rlab)
boxplot(ree,names=c(expression(mu["VA"[0]]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]^{2}),expression(omega["k"]^{2}),expression(omega[beta]^{2}),expression(omega["E"["max"]]^{2}),expression(sigma["i"])))
#bplot(ree,style="quantile",labels=c("µ_VA","µ_k","µ_be","µ_Emax","µ_ED50","w2_VA","w2_k","w2_be","w2_Emax","w2_ED50","sigma_i"))
abline(h=0,col='red')

# par(mfrow=c(1,1))
RB<-apply(ree,2,mean)
names(RB)<-vecName

RSE<-100 * SEemp / c(paramValue,omega,sigma_i)

source("D:\\yuxin\\AlgoMultip\\FIMmixte.r")
fimVA<-funFIMpop("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,168,644,672),c(0,7,644,672),c(0,168,196,672)),c(2,2,2,1,2),c(0,100,1000),c(92,59,149)/300,300)
RSEfim<-fimVA[[5]]

RSE2<-rbind(RSE,RSEfim)
barplot(RSE2,beside = T,legend=T,ylim=c(0,60),names=c(expression(mu["VA"[0]]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]^{2}),expression(omega["k"]^{2}),expression(omega[beta]^{2}),expression(omega["E"["max"]]^{2}),expression(sigma["i"])))

par(mfrow=c(2,2))
for(i in 1:10){
  hist(gg[i,],main=paste0("histogram de ",vecName[i]))
}
