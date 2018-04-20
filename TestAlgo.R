source("D:\\yuxin\\AlgoMultip\\AlgoGo.R")

w<-AlgoMult("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0.99,0,c(1),10000)
w
plot(w[[1]],ylim=c(0,1))

res<-AlgoMult("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),0.99,c(100),c(1),3000)
res

plot(res[[1]],ylim=c(0,1))

dd<-Algo
  funFIMsm("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),0,c(1))
