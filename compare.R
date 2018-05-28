

source("D:\\yuxin\\AlgoMultip\\FIMmixte.R")
source("D:\\yuxin\\AlgoMultip\\FIMmixte_1.R")

funFIMem("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),c(2,2,2,1,2),c(500),c(1),300)

funFIMem_1("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),c(2,2,2,1,2),c(500),c(1),300)

microbenchmark(funFIMem("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),c(2,2,2,1,2),c(500),c(1),300)
,funFIMem_1("VA+(1-exp(-k*t))*(emax*dose/(ED50+dose)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(c(0,7,seq(28,672,28))),c(2,2,2,1,2),c(500),c(1),300)
,times=1000)



source("D:\\yuxin\\AlgoMultip\\MultiPop.R")
source("D:\\yuxin\\AlgoMultip\\test\\partiel.R")

microbenchmark(MultiPop("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),list(),c(2,2,2),0.99,0,c(1),c(3),0.0001,5000)
,partiel("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),list(),c(2,2,2),0.99,0,c(1),c(3),0.0001,5000)
,times=100)

library(parallel)
ncores = detectCores()
l = makeCluster(ncores-1)
cluster <- makePSOCKcluster(3)


