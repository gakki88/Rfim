source("D:\\yuxin\\AlgoMultip\\FIM.R")
source("D:\\yuxin\\AlgoMultip\\FIMmixte.r")

it<-30
t<-exp(c(-50:50)/10)
parnum<-4

a<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(0,0),list(exp(c(-50:50)/10)),0,c(1),rep(1/length(t),length(t)))[[1]]											# calculate derivatives
a<-a[,-c(5,6)]  # derivatives of Part A of the Information matrix

w<-rep(1/length(t),length(t))
  for (j in 1:3000) {						
    	
    M<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(0,0),list(exp(c(-50:50)/10)),0,c(1),w)[[2]]   # calculate Information matrix for current design
    m<-solve(M)
    d<-diag(a%*%m%*%t(a))/parnum  #D-opt
    w<-w*d^0.4												# update weights. Factor 0.99 assures convergence (similar to lambda in yu2009)
    #plot(w,ylim=c(0,1))
  }



plot(w,ylim=c(0,1))
  critc<-1/t(v)%*%m%*%v											# c-criterion value
  critd<-det(M)^(1/parnum)										# D-criterion value
  return(list(w,critc,critd,m))}




