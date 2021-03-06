AlgoMult<-function(model,paramName,paramValue,sigma,design_Init,lambda,dose,PropSubjects,nbTot,delta,iteration){
  source("D:\\yuxin\\AlgoMultip\\FIM.R")
  
  w<-rep(1/length(design_Init[[1]]),length(design_Init[[1]]))
  
  
  MFi<-list()
  for(l in 1:length(design_Init[[1]])){
    MFi[[l]] <-funFIMsm(model,paramName,paramValue,sigma,design_Init[[1]][l],dose,PropSubjects,1)
    
  }
  
  it<-0
  for(k in 1:iteration){
    #Calculate the sum of matrix with weight associated
    #idea:   Mw<-Mw+w[i]*MFi[[i]]
    MFw<-Map("*",MFi,w)     #associate weights
    Mw<-Reduce("+",MFw)     #sum of all matrices
    dm<-dim(Mw)[1]
    critd<-det(Mw)^(1/dm)                   # D-criterion value
    Dphi <- critd * solve(Mw)/dm            # calculate derivatives of function phi_D
    
#     ncores = detectCores()
#     cl = makeCluster(ncores-1)
#     clusterExport(cl, "Dphi")
    d<-sapply(MFi, FUN = function(x) sum(diag(Dphi %*% x)))         #the vector of multiplier
    #stopCluster(cl)
    
    w<- w*d^lambda/sum(w*d^lambda)          #develop weight

    if(max(d)<(1+delta)*sum(w*d)){break}  #stop criterion
    
    it<-it+1
  }
  
  
  v<-which(w>mean(w))   #by default,choose w bigger than the mean value
  ff<-funFIMsm(model,paramName,paramValue,sigma,list(design_Init[[1]][v]),dose,PropSubjects,1)
  
  e1<-(funFIMsm(model,paramName,paramValue,sigma,design_Init,dose,PropSubjects,1))[[3]]
  e2<-ff[[3]]
  eff<-e1/e2
  
  cat("\n\n*******************************************************\n")
  cat("********************** Result of Algorithm ****************\n")
  cat("\n WEIGHT **********************************\n",w,
      "\n\n IMPORTANT WEIGHT INDEX ******************\n", v,
      "\n\n IMPORTANT WEIGHT VALUE ******************\n",design_Init[[1]][v],
      "\n\n NUMBER OF INTERATIONS ******************\n",it,
      "\n\n D-CRITERION ****************************\n",critd*length(v),
      "\n\n EFFICIENCY******************************\n",eff)  
  
  cat("\n\n************************************************************\n")
  cat("******************* FISHER INFORMATION MATRIX ******************\n")
  print(ff[[1]])

  cat("\n\n******************* DETERMINANT OF THE MATRIX ******************\n", ff[[2]],
      "\n\n******************* CRITERION ******************\n",ff[[3]],
      "\n\n******************* STANDARD ERROR ******************\n",ff[[4]],
      "\n\n******************* RELATIVE STANDARD ERROR ******************\n",ff[[5]])

  plot(w,ylim=c(0,1))
  return(list(w,design_Init[[1]][v],eff))
}