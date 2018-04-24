#MultiPop<-function(model,paramName,paramValue,sigma,design_Init,lambda,dose,nbSubjects,iteration){
  #source("D:\\yuxin\\AlgoMultip\\FIM.R")
  source("D:\\insa\\inserm\\FIMmixte.R")

d<-seq(0,5000,100)
Designs<-combn(d,5)
#Designs[,1]

funFIMem("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(d),2,c(0),c(1))

w<-rep(1/dim(Designs)[2],dim(Designs)[2])

MFi<-list()
for(l in 1:dim(Designs)[2]){
  MFi[[l]] <-funFIMem("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(Designs[,l]),2,c(0),c(1))[[1]]
  
}

it<-0
for(k in 1:300){
  mw<-list()
  # if((sigma[1]==0 && sigma[2]!=0) || (sigma[2]==0 && sigma[1]!=0)){
  #   Mdim=length(paramName)+1
  # }
  # if(sigma[1]==0 && sigma[2]==0){
  #   Mdim=length(paramName)
  # }
  # if(sigma[1]!=0 && sigma[2]!=0){
  #   Mdim=length(paramName)+2
  # }
  Mdim=7
  Mw<-matrix(rep(0),nrow=Mdim,ncol=Mdim)
  for(i in 1:dim(Designs)[2]){
    mw[[i]]<-w[i]*MFi[[i]]
    Mw<-Mw+mw[[i]]
  }
  #m<-solve(Mw)
  critd<-det(Mw)^(1/Mdim)     # D-criterion value
  Dphi <- critd * solve(Mw)/Mdim            # calculate derivatives of function phi_D
  d<-c()
  for(j in 1:dim(Designs)[2]){
    d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
  }
  w1<-w
  lambda=0.9
  w<- w*d^lambda/sum(w*d^lambda)
  # if((w-w1)<1e-09){break}
  # it<-it+1
}
#plot(w,ylim=c(0,1))
#}
