
#Test for individual case
source("D:\\yuxin\\AlgoMultip\\AlgoMult.R")
w<-AlgoMult("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(1,0),list(exp(c(-50:50)/10)),0.99,0,c(1),1,0.000001,10000)
w
tail(sort(w[[1]]),4)
plot(w[[1]],ylim=c(0,1))

res<-AlgoMult("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),0.99,c(100),c(1),32,0.0001,3000)
res

plot(res[[1]],ylim=c(0,1))
plot(res[[2]])

dd<-AlgoMult("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),0.99,0,c(1),1,0.0001,5000)
plot(dd[[1]],ylim=c(0,1))
dd



#Test for population case
source("D:\\yuxin\\AlgoMultip\\MultiPop.R")
b<-MultiPop("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),list(c(1,0),c(4,672)),c(2,2,2,1,2),0.99,c(0,25,50,100,150,300,500),c(1),c(4),0.0001,5000)

b<-MultiPop("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),list(),c(2,2,2,1,2),0.99,c(0,25,50,100,150,300,500),c(1),c(4),0.0001,5000)


dos<-MultiPop("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),c(2,2,2),0.99,0,c(1),c(3),0.0001,5000)
v<- which(dos[[1]]>0.1)
dos[[2]][,v]
plot(dos[[1]],ylim=c(0,1))


pk<-MultiPop("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.7,0.02,0.06),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),c(2,2,2),0.99,100,c(1),c(4),0.00001,5000)
plot(pk[[1]],ylim=c(0,1))

abc<-MultiPop("c+(de-c)/(1+exp(b*log(t)-b*log(a)))",c("c","de","b","a"),c(0,1,1,1),c(0.09,0.09,0.09),c(1,0),list(exp(c(-50:50)/10)),c(2,2,2,2),0.99,100,c(1),c(4),0.0001,500)







source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
funFIMem("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,140,672),c(0,7,140,672),c(0,7,140,672),c(0,7,140,672)),c(2,2,2,1,2),c(0,150,300,500),c(0.25,0.25,0.25,0.25),300)

funFIMem("exp(VA)+(1-exp(-exp(k)*t))*(emax*dose/(exp(ED50)+dose)-exp(be)*exp(VA))",c("VA","k","be","emax","ED50"),c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,140,672),c(0,7,140,672),c(0,7,140,672),c(0,7,140,672)),c(1,1,1,1,1),c(0,150,300,500),c(0.25,0.25,0.25,0.25),300)

MyM<-funFIMem("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),
         list(c(0,7,140,169)),c(2,2,2,1,2),c(300),c(0.25,0.25,0.25),300)

funFIMem("exp(VA)+(1-exp(-exp(k)*t))*(emax*dose/(exp(ED50)+dose)-exp(be)*exp(VA))",c("VA","k","be","emax","ED50"),c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28),0),
         list(c(0,7,140,672),c(0,7,140,672),c(0,7,140,672),c(0,7,140,672)),c(1,1,1,1,1),c(0,150,300,500),c(0.25,0.25,0.25,0.25),300)


funFIMem("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.7,0.02,0.06),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),c(2,2,2),c(100),c(1),32)

source("D:\\yuxin\\AlgoMultip\\FIM.R")
f1<-funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 2, 9, 24, 48, 96),c(1, 6, 12, 36, 72, 120)),c(100),c(0.5,0.5),32)
