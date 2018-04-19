source("D:\\insa\\inserm\\AlgoMultip\\AlgoMult.R")

w<-AlgoMult("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0.99,0,c(1),500)
w
plot(w,ylim=c(0,1))

w1<-AlgoMult("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),0.99,c(100),c(1),100)
w1
plot(w1,ylim=c(0,1))
