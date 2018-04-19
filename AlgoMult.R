AlgoMult<-function(model,paramName,paramValue,sigma,design_Init,lambda,dose,nbSubjects,iteration){
#source("D:\\yuxin\\AlgoMultip\\FIM.R")
source("D:\\insa\\inserm\\FIM.R")

#it<-1000
#t<-exp(c(-50:50)/10)
#parnum<-4
w<-rep(1/length(design_Init[[1]]),length(design_Init[[1]]))

#MF<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0,c(1))

MFi<-list()
for(l in 1:length(design_Init[[1]])){
  MFi[[l]] <-funFIMsm(model,paramName,paramValue,sigma,design_Init[[1]][l],dose,nbSubjects)
  
}

for(k in 1:iteration){
  mw<-list()
  if((sigma[1]==0 && sigma[2]!=0) || (sigma[2]==0 && sigma[1]!=0)){
    Mdim=length(paramName)+1
  }
  if(sigma[1]==0 && sigma[2]==0){
    Mdim=length(paramName)
  }
  if(sigma[1]!=0 && sigma[2]!=0){
    Mdim=length(paramName)+2
  }
  Mw<-matrix(rep(0),nrow=Mdim,ncol=Mdim)
  for(i in 1:length(design_Init[[1]])){
    mw[[i]]<-w[i]*MFi[[i]]
    Mw<-Mw+mw[[i]]
  }

  m<-solve(Mw)
  critd<-det(Mw)^(1/length(paramName))     # D-criterion value
  Dphi <- critd * solve(Mw)/length(paramName)             # calculate derivatives of function phi_D
  d<-c()
  for(j in 1:length(design_Init[[1]])){
    d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
  }

  w<- w*d^lambda/sum(w*d^lambda)
}
return(w)
#plot(w,ylim=c(0,1))

}