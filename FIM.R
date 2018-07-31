funFIMind <- function(equation,parameters,beta,sigma,t_indiv,d,PropSubjects,nbTot){
#Name of the fixed effects parameters and variable sampling times
ParametersInEquation<-c(parameters,"t")
#All parameters to be considered
PSI<-c(parameters,"sig.inter","sig.slope")
lpsi <- length(PSI)
lengthParameters<-length(parameters)

#Prepare to compute values of model equation with fixed effect parameters
equatf <- parse(text = equation, n=-1)
f<-function(ParametersInEquation){eval(equatf[[1]])}

#residual error model
Vmodel<- parse(text = paste0("( sig.inter + sig.slope * fixed )^2"))

dose<-d

#Fixed effects parameters values
for(i in 1:lengthParameters){
  assign(parameters[i],beta[i])
}

#Parameters of standard deviation model of residual error (sig.inter+sig.slope*f)^2:
sig.inter<-sigma[1]
sig.slope<-sigma[2]

M_f<-list()
M_F <- matrix(rep(0),nrow=lpsi,ncol=lpsi)
#gather all groups of protocol
for(q in 1:length(t_indiv)){
  #sample times for an individual design
  t<-t_indiv[[q]]

  
  parameterValues <- c(beta,t)
  fixed<-f(parameterValues)                  #Fixed effect values
  if(length(t)==1){
    var<-(sig.inter+sig.slope*fixed)^2       #Variance of residual error for one (measurement/sampling time)/dose
  }else{
    var<-diag((sig.inter+sig.slope*fixed)^2) #Variance of resudual error for all times/doses
  }
  

  df<-deriv(equatf[[1]],PSI)
  mdf<-attributes(eval(df))$gradient
  dv<-deriv(Vmodel[[1]],PSI)
  mdv<-attributes(eval(dv))$gradient
  for(i in 1:lengthParameters){
    mdv[,i]=0
  }
  m<- list()
  for(i in (lengthParameters+1):lpsi){
    if(length(t)==1){
      m[[i-lengthParameters]] <- mdv[i]
    }else{
      m[[i-lengthParameters]] <- diag(mdv[,i])
    }
    
  }
  
  
  M_A <- t(mdf) %*% solve(var) %*% mdf 
  M_B <- matrix(rep(0),nrow=lpsi,ncol=lpsi)
  
  for(i in 1:2){
    M_B[lengthParameters+i,c(lpsi-1,lpsi)] <- sapply(m, function(x) 1/2 * sum(diag(m[[i]] %*% solve(var) %*% x %*% solve(var) )))
  }
  M_f[[q]] <- (M_A+M_B)*PropSubjects[q]
  M_F <-M_F+M_f[[q]]
}

M_F<-M_F*nbTot

#set names for vectors
fname <-c(sapply(1:lengthParameters, function(x) paste0("u_",parameters[x])),"sig.inter","sig.slope")
rownames(M_F) <- fname
colnames(M_F) <- fname

if(0 %in% sigma){
  M_F <- M_F[-c(lengthParameters+which(c(sigma)==0)),-c(lengthParameters+which(c(sigma)==0))]
  PSI <- PSI[-c(lengthParameters+which(sigma==0))]
}

if(length(t)==1){
  return(M_F)
}  else{
  
  deterFim <- det(M_F)
  SE <- sqrt(diag(solve(M_F)))
  RSE <- 100 * SE / c(beta,sigma)[which(c(beta,sigma)!=0)]
  CritereDopt <- deterFim^(1/length(c(beta,sigma)[which(c(beta,sigma)!=0)]))

  #write the output in the console
#   cat("******************* FISHER INFORMATION MATRIX ******************\n")
#   print(M_F)
#   
#   cat("\n\n******************* DETERMINANT OF THE MATRIX ******************\n", deterFim,
#       "\n\n******************* CRITERION ******************\n",CritereDopt,
#       "\n\n******************* STANDARD ERROR ******************\n",SE,
#       "\n\n******************* RELATIVE STANDARD ERROR ******************\n",RSE,"\n")
  
  #suppose just one group of design, plot the function and result points
#   t<-t_indiv[[1]]
#   x1<-c(beta,t)
#   fpoints<-f(x1)
#   plot(t,fpoints,xlab="Dose",ylab="Effect") #first plot the result points
#   t<-seq(min(t)-0.5,max(t),1)
#   x<-c(beta,t)
#   fvalue<-f(x)
#   lines(t,fvalue) #then plot the line of function
  
  return(list(M_F,deterFim,CritereDopt,SE,RSE))
}



}


#exercice1
# funFIMind("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),c(100),c(1),32)

#exercice2
# funFIMind("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.6,0.07),list(c(0.5, 2, 9, 24, 48, 96),c(1, 6, 12, 36, 72, 120)),c(100),c(0.5,0.5),32)

#emample dose-response
# funFIMind("Emax*t/(t+D50)+E0",c("Emax","D50","E0"),c(30,500,5),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),0,c(1),1)
 