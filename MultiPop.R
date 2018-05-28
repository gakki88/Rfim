http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/

MultiPop<-function(model,paramName,paramValue,omega,sigma,S_times,fixed_times,Trand,lambda,dose,PropSubjects,nbTimes,delta,iteration){
    
    source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
    #for(n in 1:length(nbTimes)){

      start_time <- Sys.time()
      Designs_ini<-combn(S_times[[1]],nbTimes[1])
      
      if(length(fixed_times)>0){
        for(f in 1:length(fixed_times)){
          Designs_ini<-Designs_ini[,which(Designs_ini[fixed_times[[f]][1],]==fixed_times[[f]][2])]
        }
      }

      MFi<-list()
      Designs<-matrix(nrow=nbTimes[1])
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
      Designs<-Designs[,-1]
      o<-o-1
      w<-rep(1/(o),(o))
      it<-0
      
      for(k in 1:iteration){
        
        #Calculate the sum of matrix with weight associated
        #idea:   Mw<-Mw+w[i]*MFi[[i]]
        MFw<-Map("*",MFi,w)     #associate weights
        Mw<-Reduce("+",MFw)     #sum of all matrices
        dm<-dim(Mw)[1]
        critd<-det(Mw)^(1/dm)                   # D-criterion value
        Dphi <- critd * solve(Mw)/dm            # calculate derivatives of function phi_D
        
        ncores = detectCores()
        cl = makeCluster(ncores-1)
        clusterExport(cl, "Dphi")
        d<-unlist(parLapply(cl,MFi, FUN = function(x) sum(diag(Dphi %*% x))))         #the vector of multiplier
        stopCluster(cl)
        
        w<- w*d^lambda/sum(w*d^lambda)          #develop weight
        
        if(max(d)<(1+delta)*sum(w*d)){break}    #stop criterion
        it<-it+1                                #number of iteration calculator
        if(it%%100==0) print(it)
      }
      end_time <- Sys.time()
      plot(w,ylim=c(0,1))
      v<-which(w>mean(w))
      w1<-w[v]
      ddose<-c(dose,dose)
      OptDoses<-ddose[v%%length(dose)+length(dose)]
      #   e1<-funFIMem(model,paramName,paramValue,omega,sigma,S_times,Trand,dose,PropSubjects,1)[[3]]
      #   e2<-c()
      #   for(q in 1:length(v)){
      #     e2<-c(e2,funFIMem(model,paramName,paramValue,omega,sigma,list(Designs[,v[q]]),Trand,dose,PropSubjects,1)[[3]])
      #   }
      #   
      #   cat( "\n WEIGHT***********************************\n",w)
      cat("\n NUMBER OF ITERATION**********************\n",it,
          "\n UPPER INDEX*****************************\n",v,
          "\n UPPER WEIGHTS***********************\n",w1,
          "\n DOSES CORRESPONDANT******************\n",OptDoses)
      #       "\n DETERMINANT**************************\n",e2,
      #       "\n EFFICIENCY**************************\n",e1/e2,
      cat("\n VALUED DESIGNS***********************\n")
      print(Designs[,v])
      cat("\n EXECUTE TIME********************",end_time - start_time,"\n")
      
     #}
     #return(list(w,Designs,it))
  }