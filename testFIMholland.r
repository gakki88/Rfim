#funFIMsm <- function(equation,parameter,sigma,t_indiv,d,nbSubjects){
#Name of the fixed effects parameters
#paramName<-c("ka","V","Cl")
#paramName<-c("Emax","C50","S0") 
paramName<-c("c","de","b","e")

paramF<-c(paramName,"t")
#"c+(de-c)/(1+exp(b*log(t)-b*log(e)))",c(0,1,1,1),c(0,0),list(exp(c(-50:50)/10)),0,c(1)
#model equation
equatf <- parse(text = "c+(de-c)/(1+exp(b*log(t)-b*log(e)))", n=-1)
f<-function(paramF){eval(equatf[[1]])}

c<-0  												# parameter settings
de<-1
b<-1
e<-1



#sample times for an individual design
#t <- t_indiv
M_f<-list()
M_F <- matrix(rep(0),nrow=length(PSI),ncol=length(PSI))
#gather all groups of protocol

  t<-exp(c(-50:50)/10)



beta<-c(0,1,1,1)

param <- c(beta,t)
fixed<-f(param)




df<-deriv(equatf[[1]],paramName)
mdf<-attributes(eval(df))$gradient

for(i in 1:length(paramName)){
  mdv[,i]=0
}
m<- list()
for(i in (length(paramName)+1):lpsi){
  m[[i-length(paramName)]] <- diag(mdv[,i])
}

w<-rep(1/101,101)
M_A <- t(mdf) %*% diag(w) %*% mdf 

}
# if(sig.slope ==0){
#   M_F <- M_F[,-c(lpsi)]
#   M_F <- M_F[-c(lpsi),]
#   PSI <- PSI[-c(lpsi)]
# }else if(sig.inter == 0){
#   M_F <- M_F[,-c(lpsi-1)]
#   M_F <- M_F[-c(lpsi-1),]
#   PSI <- PSI[-c(lpsi-1)]
# }


deterFim <- det(M_F)
# SE <- diag(sqrt(solve(M_F)))
# RSE <- 100 * SE / c(beta,sigma)
CritereDopt <- deterFim^(1/length(PSI))

#set names for vectors -- time protocol
# rownames(M_F) <- c("ka_fixed","V_fixed","Cl_fixed","sig.inter","sig.slope")
# colnames(M_F) <- c("ka_fixed","V_fixed","Cl_fixed","sig.inter","sig.slope")
# names(SE) <- c("ka_fixed","V_fixed","Cl_fixed","sig.inter","sig.slope")
# names(RSE) <- c("ka_fixed","V_fixed","Cl_fixed","sig.inter","sig.slope")

#set names -- dose response
# rownames(M_F) <- c("Emax_fixed","C50_fixed","S0_fixed","sig.inter","sig.slope")
# colnames(M_F) <- c("Emax_fixed","C50_fixed","S0_fixed","sig.inter","sig.slope")
# names(SE) <- c("Emax_fixed","C50_fixed","S0_fixed","sig.inter","sig.slope")
# names(RSE) <- c("Emax_fixed","C50_fixed","S0_fixed","sig.inter","sig.slope")

#return(list(M_F,deterFim,SE,RSE,CritereDopt))
#write the output into a text file
#sink('D:\\insa\\inserm\\ExoPFIM\\ex2\\rfimI.txt')
cat("******************* FISHER INFORMATION MATRIX ******************\n")
print(M_F)

cat("\n\n******************* DETERMINANT OF THE MATRIX ******************\n", deterFim,
    "\n\n******************* STANDARD ERROR ******************\n",SE,
    "\n\n******************* RELATIVE STANDARD ERROR ******************\n",RSE,
    "\n\n******************* CRITERION ******************\n",CritereDopt)

#sink()

}
#exercice1
#funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c(1.6,8,0.13),c(0.6,0.07),c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120),c(100),c(32))

#exercice2
#funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 2, 9, 24, 48, 96),c(1, 6, 12, 36, 72, 120)),c(100),c(16,16))

#emample dose-response
#funFIMsm("Emax*t/(t+C50)+S0",c(30,500,5),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),0,c(1))

#FIM de Holland
funFIMsm("c+(de-c)/(1+exp(b*log(t)-b*log(e)))",c(0,1,1,1),c(0,0),list(exp(c(-50:50)/10)),0,c(1))

