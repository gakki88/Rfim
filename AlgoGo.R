source("D:\\yuxin\\AlgoMultip\\FIM.R")

it<-1000
t<-exp(c(-50:50)/10)
parnum<-4
w<-rep(1/length(t),length(t))

#MF<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0,c(1))

MFi<-list()
for(l in 1:length(t)){
  MFi[[l]] <-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(t[l]),0,c(1))

}

for(k in 1:it){
  mw<-list()
  Mw<-matrix(rep(0),nrow=parnum+1,ncol=parnum+1)
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
}
w
plot(w,ylim=c(0,1))