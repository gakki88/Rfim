MultiPop("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0, 7, seq(28,672,28))),c(2,2,2,1,2),0.99,c(0, 25, 50, 100, 150, 300, 500 ),c(1),c(4),0.0001,500)


MultiPop<-function("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0, 7, seq(28,672,28))),c(2,2,2,1,2),0.99,c(0, 25, 50, 100, 150, 300, 500 ),c(1),c(4),0.0001,500){
  source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
  
  
  
    start_time <- Sys.time()
    Designs_ini<-combn(c(0, 7, seq(28,672,28)),4)
    dose<-c(0, 25, 50, 100, 150, 300, 500 )
    MFi<-list()
#     MFi2<-list()
#     MFi3<-list()
#     MFi4<-list()
#     MFi5<-list()
#     MFi6<-list()
#     MFi7<-list()
    Designs<-matrix(nrow=4)
    o<-1
    print(dim(Designs_ini)[2])
    for(l in 1:dim(Designs_ini)[2]){
      #tryCatch({
      for(m in 1:length(dose)){
      #the number of design groups need to change in accordance with the number of doses
      MFi[[o]] <-funFIMem("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(Designs_ini[,l]),c(2,2,2,1,2),dose[m],c(1),1)[[1]]
      Designs<-cbind(Designs,Designs_ini[,l])
      o<-o+1
      
#       MFi2[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[2],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
#       
#       MFi3[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[3],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
#       
#       MFi4[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[4],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
#       
#       MFi5[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[5],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
#       
#       MFi6[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[6],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
#       
#       MFi7[[o]] <-funFIMem(model,paramName,paramValue,omega,sigma,list(Designs_ini[,l]),Trand,dose[7],PropSubjects,1)[[1]]
#       Designs<-cbind(Designs,Designs_ini[,l])
#       o<-o+1
      }
      
      #}, error=function(e){})
      if(l%%100==0) print(l)
    }
    Designs<-Designs[,-1]
    #dd<-dim(Designs)[2]*length(dose)
    #w<-rep(1/dd,dd)
    o<-o-1
    w<-rep(1/(o),(o))
    it<-0
    print(o)
    for(k in 1:50){
      
      #Initialize the matrix with the dimension above
      dm<-dim(MFi[[1]])[1]
      Mw<-matrix(rep(0),nrow=dm,ncol=dm)
      for(i in 1:o){
        Mw<-Mw+w[i]*MFi[[i]]
#         Mw<-Mw+w[i]*MFi2[[i]]
#         Mw<-Mw+w[i]*MFi3[[i]]
#         Mw<-Mw+w[i]*MFi4[[i]]
#         Mw<-Mw+w[i]*MFi5[[i]]
#         Mw<-Mw+w[i]*MFi6[[i]]
#         Mw<-Mw+w[i]*MFi7[[i]]
      }
      critd<-det(Mw)^(1/dm)     # D-criterion value
      Dphi <- critd * solve(Mw)/dm            # calculate derivatives of function phi_D
      d<-c()
      for(j in 1:o){
        d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi2[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi3[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi4[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi5[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi6[[j]])))
#         d<-c(d,sum(diag(Dphi %*% MFi7[[j]])))
      }
      lambda<-0.99
      w<- w*d^lambda/sum(w*d^lambda)
      delta<-0.0001
      if(max(d)<(1+delta)*sum(w*d)){break}
      it<-it+1
      print(it)
    }
    end_time <- Sys.time()
    plot(w,ylim=c(0,1))
    v<-which(w>mean(w)*600)
    w1<-w[v]
    OptDoses<-dose[v%%length(dose)+1]
  
    #   cat( "\n WEIGHT***********************************\n",w)
    cat("\n NUMBER OF ITERATION**********************\n",it,
        "\n UPPER INDEX*****************************\n",v,
        "\n UPPER WEIGHTS***********************\n",w1,
        "\n DOSES CORRESPONDANT******************\n",OptDoses)

    cat("\n VALUED DESIGNS***********************\n")
    print(Designs[,v])
    cat("\n EXECUTE TIME********************",end_time - start_time,"\n")
    
  
  #return(list(w,Designs,it))
}
