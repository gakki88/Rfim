#Name of the fixed effects parameters
paramName<-c("ka","V","Cl")

paramF<-c(paramName,"t")

#model equation
form1<- c("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))")


form2<-paste0("( sig.inter + sig.slope * ", form1,")^2")
#form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
Vmodel<- parse(text = form2)

PSI<-c(paramName,"sig.inter", "sig.slope")


#--------------------------------------------
#calculate Fisher Information Matrix
#--------------------------------------------

#number of subjects
n=32
#enter the number of groups
cat(" Enter the number of subject groups\n")
nn <- as.integer(readline(prompt = ""))

#dose value
dose<-c(100)

#enter designs
for(i in 1:nn){
  cat("Enter the ",i,"th group initial design \n")
  cat("Enter the quantity of this design \n")
  ni <- as.integer(readline(prompt = ""))
  ti=c()
  for(j in 1:ni){
    tij <- as.numeric(readline(prompt = ""))
    ti <- c(ti,tij)
  }
  
  
  
}
t<-c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)

#Fixed effects parameters values
beta<-c(1.6,8,0.13)

#(Diagonal Matrix of) variance for inter-subject random effects:
omega<-diag(c(0.7,0.02,0.06))


#Random effect model (1) = additive  (2) = exponential 
#------------------------------------------------------------------
Trand<-2;

if ( Trand == 1 ) {
  beta <- beta + rnorm(length(t)*length(paramName), mean = 0, sd = omega)
  bm <- matrix(beta,nrow=3,ncol=12)
  for(i in 1:length(paramName)){
    form11 <- gsub(paramName[1], paste0("(",paramName[1],"+b)"), form1)
    
  }
  
}
else if ( Trand == 2){
  b = rnorm(length(t)*length(paramName), mean = 0, sd = omega)
  beta <- beta * t(exp(b))
  bm <- matrix(beta,nrow=3,ncol=12)
  for(i in 1:length(paramName)){
    form11 <- gsub(paramName[1], paste0("(",paramName[1],"*exp(b))"), form1)
  }
}


#calculate matrix E for n individuals
equatf <- parse(text = form11, n=-1)
f<-function(paramF){eval(equatf[[1]])}

ka = beta[1]
V = beta[2]
Cl = beta[3]
param <- list()
for(i in 1:nn){
  param <- c(param,beta,t[i])
}

fixed<-f(param)



#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
sig.inter<-0.6
sig.slope<-0.07


#lambda<-c(sig.inter, sig.slope)

#psi<-c(beta,lambda)

var<-diag((sig.inter+sig.slope*fixed)^2)

epsilon<-rnorm(length(t), mean = 0, sd = diag(var))

obs<-fixed+epsilon


df<-deriv(equatf[[1]],PSI)
mdf<-attributes(eval(df))$gradient
dv<-deriv(Vmodel[[1]],PSI)
mdv<-attributes(eval(dv))$gradient
for(i in 1:length(paramName)){
  mdv[,i]=0
}
m<- list()
for(i in (length(paramName)+1):length(PSI)){
  m[[i-length(paramName)]] <- diag(mdv[,i])
}


M_A <- t(mdf) %*% solve(var) %*% mdf 
M_B <- matrix(rep(0),nrow=length(PSI),ncol=length(PSI))
for(i in 1:2){
  for(j in 1:2){
    M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag((m[[i]] %*% solve(var) %*% m[[j]] %*% solve(var))))
  }
}
M_F <- M_A+M_B
