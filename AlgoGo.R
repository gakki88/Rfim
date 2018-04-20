AlgoMult<-function(model,paramName,paramValue,sigma,design_Init,lambda,dose,nbSubjects,iteration){
  source("D:\\yuxin\\AlgoMultip\\FIM.R")
  cat("********************Avant optimiser*******************\n")
  funFIMsm(model,paramName,paramValue,sigma,design_Init,dose,nbSubjects)
  w<-rep(1/length(design_Init[[1]]),length(design_Init[[1]]))
  
  
  MFi<-list()
  for(l in 1:length(design_Init[[1]])){
    MFi[[l]] <-funFIMsm(model,paramName,paramValue,sigma,design_Init[[1]][l],dose,nbSubjects)
    
  }
  
  it<-0
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
    #m<-solve(Mw)
    critd<-det(Mw)^(1/Mdim)     # D-criterion value
    Dphi <- critd * solve(Mw)/Mdim            # calculate derivatives of function phi_D
    d<-c()
    for(j in 1:length(design_Init[[1]])){
      d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
    }
    w1<-w
    w<- w*d^lambda/sum(w*d^lambda)
    if((w-w1)<1e-09){break}
    it<-it+1
  }
  
  v<-which(w>1/length(paramName)*0.7)
  cat("\n\n*********************Apres optimiser****************\n")
  funFIMsm(model,paramName,paramValue,sigma,list(design_Init[[1]][v]),dose,nbSubjects)
  
  cat("\n**************** Result of Algorithm ****************\n")
  
  #return(list(w,v,design_Init[[1]][v],it))
  cat("******************* WEIGHT ******************\n",w,
      "\n\n******************* IMPORTANT WEIGHT INDEX ******************\n", v,
      "\n\n******************* IMPORTANT WEIGHT VALUE ******************\n",design_Init[[1]][v],
      "\n\n******************* NUMBER OF INTERATIONS ******************\n",it,
      "\n\n******************* D-CRITERION ******************\n",critd*length(v))
}
