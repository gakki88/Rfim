MultiPop<-function(model,paramName,paramValue,omega,sigma,design_Init,Trand,lambda,dose,PropSubjects,nbD,delta,iteration){
  source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
  
  
  
for(n in 1:length(nbD)){
  start_time <- Sys.time()
  Designs_ini<-combn(design_Init[[1]],nbD[n])
  
  MFi<-list()
  Designs<-matrix(nrow=nbD[n])
  for(l in 1:dim(Designs_ini)[2]){
    tryCatch({
    MFi[[l]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose,PropSubjects,1)[[1]]
    Designs<-cbind(Designs,Designs_ini[,l])
    }, error=function(e){})
    
  }
  Designs<-Designs[,-1]
  dd<-dim(Designs)[2]
  w<-rep(1/dd,dd)
  it<-0
  
  #Prepare the dimension for the matrix associated with weights
  if((sigma[1]==0 && sigma[2]!=0) || (sigma[2]==0 && sigma[1]!=0)){
    Mdim=length(paramName)*2+1
  }
  if(sigma[1]!=0 && sigma[2]!=0){
    Mdim=length(paramName)*2+2
  }
  for(k in 1:iteration){

    #Initialize the matrix with the dimension above
    Mw<-matrix(rep(0),nrow=Mdim,ncol=Mdim)
    
    for(i in 1:dd){
      Mw<-Mw+w[i]*MFi[[i]]
    }
    
    critd<-det(Mw)^(1/Mdim)     # D-criterion value
    Dphi <- critd * solve(Mw)/Mdim            # calculate derivatives of function phi_D
    d<-c()
    for(j in 1:dd){
      d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
    }
    w<- w*d^lambda/sum(w*d^lambda)
    if(max(d)<(1+delta)*sum(w*d)){break}
    it<-it+1
  }
  end_time <- Sys.time()
  plot(w,ylim=c(0,1))
  v<-which(w>mean(w))
  w1<-w[v]
  e1<-funFIMem(model,paramName,paramValue,omega,sigma,design_Init,Trand,dose,PropSubjects,1)[[3]]
  e2<-c()
  for(q in 1:length(v)){
    e2<-c(e2,funFIMem(model,paramName,paramValue,omega,sigma,list(Designs[,v[q]]),Trand,dose,PropSubjects,1)[[3]])
  }
  
  # "\n WEIGHT***********************************\n",w,
  cat("\n NUMBER OF ITERATION**********************\n",it,
      "\n UPPER INDEX*****************************\n",v,
      "\n UPPER WEIGHTS***********************\n",w1,
      "\n DETERMINANT**************************\n",e2,
      "\n EFFICIENCY**************************\n",e1/e2,
      "\n VALUED DESIGNS***********************\n")
  print(Designs[,v])
  cat("\n EXECUTE TIME********************",end_time - start_time,"\n")

}
  #return(list(w,Designs,it))
}
