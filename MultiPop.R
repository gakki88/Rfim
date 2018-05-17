

MultiPop<-function(model,paramName,paramValue,omega,sigma,design_Init,Trand,lambda,dose,PropSubjects,nbTimes,delta,iteration){
  source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
  
  
  
for(n in 1:length(nbTimes)){
  start_time <- Sys.time()
  Designs_ini<-combn(design_Init[[1]],nbTimes[n])
  
  MFi<-list()
  Designs<-matrix(nrow=nbTimes[n])
  o<-1
  for(l in 1:dim(Designs_ini)[2]){
    #tryCatch({
    for(m in 1:length(dose)){
      #the number of design groups need to change in accordance with the number of doses
      MFi[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[m],PropSubjects,1)[[1]]
      Designs<-cbind(Designs,Designs_ini[,l])
      o<-o+1
    }
      
    #}, error=function(e){})
    
  }
  Designs<-Designs[,-1]
  #dd<-dim(Designs)[2]*length(dose)
  #w<-rep(1/dd,dd)
  o<-o-1
  w<-rep(1/(o),(o))
  it<-0
  
  for(k in 1:iteration){

    #Initialize the matrix with the dimension above
    dm<-dim(MFi[[1]])[1]
    Mw<-matrix(rep(0),nrow=dm,ncol=dm)
    for(i in 1:o){
      Mw<-Mw+w[i]*MFi[[i]]
    }
    critd<-det(Mw)^(1/dm)     # D-criterion value
    Dphi <- critd * solve(Mw)/dm            # calculate derivatives of function phi_D
    d<-c()
    for(j in 1:o){
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
#   e1<-funFIMem(model,paramName,paramValue,omega,sigma,design_Init,Trand,dose,PropSubjects,1)[[3]]
#   e2<-c()
#   for(q in 1:length(v)){
#     e2<-c(e2,funFIMem(model,paramName,paramValue,omega,sigma,list(Designs[,v[q]]),Trand,dose,PropSubjects,1)[[3]])
#   }
#   
#   cat( "\n WEIGHT***********************************\n",w)
   cat("\n NUMBER OF ITERATION**********************\n",it,
       "\n UPPER INDEX*****************************\n",v,
      "\n UPPER WEIGHTS***********************\n",w1)
#       "\n DETERMINANT**************************\n",e2,
#       "\n EFFICIENCY**************************\n",e1/e2,
      cat("\n VALUED DESIGNS***********************\n")
  print(Designs[,v])
   cat("\n EXECUTE TIME********************",end_time - start_time,"\n")

}
  #return(list(w,Designs,it))
}
