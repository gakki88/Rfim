# source("D:\\yuxin\\AlgoMultip\\FIMam.R")
source("D:\\yuxin\\code\\FIM.R")

it<-300
t<-c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)
parnum<-3
w<-rep(1/length(t),length(t))

a<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[1]]  				# calculate derivatives
a<-a[,-c(parnum+1,parnum+2)]  # derivatives of Part A of the Information matrix

MF<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[2]]   # calculate Information matrix for current design
vv<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[3]]
mm<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[4]]
mv<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[5]]
mb<- funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(t),c(100),c(1))[[6]]

w<-rep(1/length(t),length(t))
for(k in 1:100){
  M_A <- matrix(rep(0),nrow=parnum+2,ncol=parnum+2)
  M_A[1:parnum,1:parnum] <-  t(a) %*% solve(vv) %*% diag(w) %*% a 
  M_B <- matrix(rep(0),nrow=parnum+2,ncol=parnum+2)
  for(i in 1:2){
    for(j in 1:2){
      M_B[parnum+i,parnum+j] <- 1/2 * sum(diag((mm[[i]]%*% diag(w) %*% solve(vv)%*%mm[[j]] %*% solve(vv))))
    }
  }
  Mw <-M_A +M_B
  
  m<-solve(Mw)
  critd<-det(Mw)^(1/parnum)     # D-criterion value
  Dphi <- critd * solve(Mw)/parnum             # calculate derivatives of function phi_D
  
  d1<-a%*%Dphi[1:parnum,1:parnum]%*%t(a)%*%solve(vv)
  dd<-diag((mm[[1]] %*% solve(vv)%*%mm[[1]] %*% solve(vv)))*Dphi[parnum+1,parnum+1]+diag((mm[[1]] %*% solve(vv)%*%mm[[2]] %*% solve(vv)))*Dphi[parnum+2,parnum+1]+diag((mm[[2]] %*% solve(vv)%*%mm[[1]] %*% solve(vv)))*Dphi[parnum+1,parnum+2]+diag((mm[[2]] %*% solve(vv)%*%mm[[2]] %*% solve(vv)))*Dphi[parnum+2,parnum+2]
  d<-diag(d1)+dd/96
  #w<- w*diag(d1)^0.99/sum(w*diag(d1))
  w<- w*d^0.99/sum(w*d)
}
w
plot(w,ylim=c(0,1))

w1<-rep(1/length(t),length(t))   
for (j in 1:3000) {    				
  
  #   M<-funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(0,0),list(exp(c(-50:50)/10)),0,c(1))[[2]]   # calculate Information matrix for current design
  M<-t(a)%*%diag(w1)%*%a
  mh<-solve(M)
  d<-diag(a%*%mh%*%t(a))/parnum  #D-opt
  w1<-w1*d^0.99												# update weights. Factor 0.99 assures convergence (similar to lambda in yu2009)
  
}
plot(w1,ylim=c(0,1))
