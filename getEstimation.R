getEstimation<-function(N){
  
  directory <-paste0("D:\\yuxin\\Simulation3\\resultCDC4\\md_",1,"\\project0407")
#   directory <-paste0("C:\\yuxinCDC\\resultCDC3\\md_",1,"\\project3")
  setwd(directory)
  Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
  datas<-Data[,2]
  
  for(i in 2:N){
    tryCatch({
      directory <-paste0("D:\\yuxin\\Simulation3\\resultCDC4\\md_",i,"\\project0407")
#       directory <-paste0("C:\\yuxinCDC\\resultCDC3\\md_",i,"\\project3")
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
vecName<-c("�_VA0","�_k","�_be","�_Emax","�_ED50","w2_VA","w2_k","w2_be","w2_Emax","sigma_i")
#vecName<-c(expression(mu["VA"]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]^{2}),expression(omega["k"]^{2}),expression(omega[beta]^{2}),expression(omega["E"["max"]]^{2}),expression(sigma["i"]))
SEemp<-c()
for(i in 1:5){
  SEemp <- c(SEemp,sd(gg[i,]))
  # boxplot(gg[i,])
}
for(i in 6:9){
  SEemp <- c(SEemp,sd(sqrt(gg[i,])))
  # boxplot(gg[i,])
}
SEemp <- c(SEemp,sd(gg[10,]))

names(SEemp) <- vecName

paramValue<-c(55,0.005,0.2,30,150)
omega<-sqrt(c(0.07,0.5,1,150))
sigma_i<-sqrt(28)
ree<-c((gg[1,]-paramValue[1])/paramValue[1]*100)
ree2<-ree
for(i in 2:5){
  ree<-c(ree,(gg[i,]-paramValue[i])/paramValue[i]*100)
  ree2<-cbind(ree2,(gg[i,]-paramValue[i])/paramValue[i]*100)
}
for(i in 1:4){
  ree<-c(ree,(sqrt(gg[5+i,])-omega[i])/omega[i]*100)
  ree2<-cbind(ree2,(sqrt(gg[5+i,])-omega[i])/omega[i]*100)
}
ree<-c(ree,(gg[10,]-sigma_i)/sigma_i*100)
ree2<-cbind(ree2,(gg[10,]-sigma_i)/sigma_i*100)
ree1 <- data.frame(values=ree, parameters = rep(vecName, each = 500))
#colnames(ree)<-vecName
#Boxplots ree#######################
#library(Rlab)
# boxplot(ree2,names=c(expression(mu["VA"[0]]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]^{2}),expression(omega["k"]^{2}),expression(omega[beta]^{2}),expression(omega["E"["max"]]^{2}),expression(sigma["i"])))
#bplot(ree,style="quantile",labels=c("�_VA","�_k","�_be","�_Emax","�_ED50","w2_VA","w2_k","w2_be","w2_Emax","w2_ED50","sigma_i"))


f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
customerr <- function(x) {
  r <- quantile(x, probs = c(0.05,  0.95))
  names(r) <- c("ymin", "ymax")
  r
}

library(ggplot2) #,fill=parameters

ree1$parameters <- factor(ree1$parameters, levels = ree1$parameters)
ggplot(ree1) +aes(x=parameters, y=values,color="") +
  theme_bw(15)+
  labs(y="REE")+
  stat_summary(fun.data = f, geom='boxplot',position="dodge")+
  geom_hline( yintercept=0)+
  stat_summary(fun.data = customerr, geom="errorbar",position=position_dodge(0.9),size=1, width=0.6) +
  scale_x_discrete( labels=c(expression(mu["VA"[0]]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]),expression(omega["k"]),expression(omega[beta]),expression(omega["E"["max"]]),expression(sigma)))

#ggsave("D:\\yuxin\\Simulation2\\resultatSimulation\\ree11111.png")

#ggplot(ree1)+aes(x=parameters, y=values,fill=parameters)+  geom_boxplot(fill='#A4A4A4', color="darkred")
#ggsave("D:\\yuxin\\Simulation2\\resultatSimulation\\ree11.png")
#abline(h=0,col='red')
#ggsave("D:\\yuxin\\Simulation2\\resultatSimulation\\ree1.png")

# par(mfrow=c(1,1))
ree2<-ree[1:500]
for(i in 2:10)
{
  ree2<-cbind(ree2,ree[(500*(i-1)+1):(500*i)])
}
RB<-apply(ree2,2,mean)
names(RB)<-vecName

RSE<-100 * SEemp / c(paramValue,omega,sigma_i)

source("D:\\yuxin\\AlgoMultip\\FIMmixte.r")
fimVA<-funFIMpop("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,168,196,672),c(0,7,644,672),c(0,7,644,672),c(0,168,196,672)),c(2,2,2,1,2),c(0,50,100,500),c(107,13,48,132)/300,300)
RSEfim<-fimVA[[5]]
RSEfim[c(6:9)]<-RSEfim[c(6:9)]/2
RSE_CTS<-RSE
RSE_FIM<-RSEfim
RSE2<-rbind(RSE_CTS,RSE_FIM)
barplot(RSE2,beside = T,col=c("darkblue","red"),legend=T,ylim=c(0,50),names=c(expression(mu["VA"[0]]),expression(mu["k"]),expression(mu[beta]),expression(mu["E"["max"]]),expression(mu["ED"[50]]),expression(omega["VA"]),expression(omega["k"]),expression(omega[beta]),expression(omega["E"["max"]]),expression(sigma)))

# par(mfrow=c(2,2))
# for(i in 1:10){
#   hist(gg[i,],main=paste0("histogram de ",vecName[i]))
# }
