source("D:\\insa\\inserm\\FIMmixte.R")
Designs_ini<-combn(c(0, 7, seq(28,672,28)),4)
dose<-c(0, 25, 50, 100, 150, 300, 500 )
MFi<-list()
Designs<-matrix(nrow=4)
o<-1
print(dim(Designs_ini)[2])
MF2<-matrix(nrow=10,ncol=10)
for(l in 1:dim(Designs_ini)[2]){
  #tryCatch({
  for(m in 1:length(dose)){
    #the number of design groups need to change in accordance with the number of doses
    MFi[[o]] <-funFIMem("VA+(1-exp(-k*t))*(emax*d/(ED50+d)-be*VA)",c("VA","k","be","emax","ED50"),c(55,0.005,0.2,30,150),c(0.07,0.5,1,150,0),c(sqrt(28),0),list(Designs_ini[,l]),c(2,2,2,1,2),dose[m],c(1),1)[[1]]
    #Designs<-cbind(Designs,Designs_ini[,l])
    o<-o+1

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
for(k in 1:500){
  
  #Initialize the matrix with the dimension above
  #dm<-dim(MFi[[1]])[1]
  #Mw<-matrix(rep(0),nrow=dm,ncol=dm)
  # for(i in 1:o){
  #   Mw<-Mw+w[i]*MFi[[i]]
  # }
  MFw<-Map("*",MFi,w)
  Mw<-Reduce("+",MFw)
  critd<-det(Mw)^(1/dm)     # D-criterion value
  Dphi <- critd * solve(Mw)/dm            # calculate derivatives of function phi_D
  # d<-c()
  # for(j in 1:o){
  #   d<-c(d,sum(diag(Dphi %*% MFi[[j]])))
  # }
  #d<-lapply(MFi, FUN = function(x) Dphi %*% x)
  d<-unlist(lapply(MFi, FUN = function(x) sum(diag(Dphi %*% x))))
  lambda<-0.99
  w<- w*d^lambda/sum(w*d^lambda)
  delta<-0.0001
  if(max(d)<(1+delta)*sum(w*d)){break}
  it<-it+1
  print(it)
}
end_time <- Sys.time()
plot(w,ylim=c(0,1))
v<-which(w>0.03)
w1<-w[v]
ddose<-c(dose,dose)
OptDoses<-ddose[v%%length(dose)+length(dose)]

#   cat( "\n WEIGHT***********************************\n",w)
cat("\n NUMBER OF ITERATION**********************\n",it,
    "\n UPPER INDEX*****************************\n",v,
    "\n UPPER WEIGHTS***********************\n",w1,
    "\n DOSES CORRESPONDANT******************\n",OptDoses)

cat("\n VALUED DESIGNS***********************\n")
print(Designs[,v])
cat("\n EXECUTE TIME********************",end_time - start_time,"\n")


#return(list(w,Designs,it))
