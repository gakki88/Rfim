

partiel<-function(model,paramName,paramValue,omega,sigma,S_times,fixed_times,Trand,lambda,dose,PropSubjects,nbTimes,delta,iteration){
  
  source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
  for(n in 1:length(nbTimes)){
    start_time <- Sys.time()
    Designs_ini<-combn(S_times[[1]],nbTimes[n])
    
    if(length(fixed_times)>0){
      for(f in 1:length(fixed_times)){
        Designs_ini<-Designs_ini[,which(Designs_ini[fixed_times[[f]][1],]==fixed_times[[f]][2])]
      }
    }
    
    MFi<-list()
    Designs<-matrix(nrow=nbTimes[n])
    o<-1
    for(l in 1:dim(Designs_ini)[2]){
      for(m in 1:length(dose)){
        #the number of design groups need to change in accordance with the number of doses
        MFi[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[m],PropSubjects,1)[[1]]
        Designs<-cbind(Designs,Designs_ini[,l])
        o<-o+1
      }
      
      if(l%%100==0) print(l)
    }
    end_time <- Sys.time()
}
}