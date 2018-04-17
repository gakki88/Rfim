source("D:\\yuxin\\AlgoMultip\\FIM.R")

it<-500
t<-c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)
parnum<-3
w<-rep(1/length(t),length(t))

#MF<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0,c(1))

MFi<-list()
for(ll in 1:length(t)){
  MFi[[ll]] <-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(ll),c(100),c(1))
  
}


for(k in 1:it){
  mw<-list()
  Mw<-matrix(rep(0),nrow=parnum+2,ncol=parnum+2)
  for(i in 1:length(t)){
    mw[[i]]<-w[i]*MFi[[i]]
    Mw<-Mw+mw[[i]]
  }
  
  m<-solve(Mw)
  critd<-det(Mw)^(1/parnum)     # D-criterion value
  Dphi <- critd * solve(Mw)/parnum             # calculate derivatives of function phi_D
  d<-c()
  for(j in 1:length(t)){
    d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
  }
  
  w<- w*d^0.99/sum(w*d^0.99)
  plot(w,ylim=c(0,1))
}
w
plot(w,ylim=c(0,1))