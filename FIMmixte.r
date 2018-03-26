#Name of the fixed effects parameters
paramName<-c("ka","V","Cl")

paramF<-c(paramName,"t")

#model equation
form1<- c("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))")


PSI<-c(paramName,"sig.inter", "sig.slope")


#--------------------------------------------
#calculate Fisher Information Matrix
#--------------------------------------------

#number of subjects
n=32
#enter the number of groups
# cat(" Enter the number of subject groups\n")
# nn <- as.integer(readline(prompt = ""))

#dose value
dose<-c(100)

#enter designs
# t <- c()
# for(i in 1:nn){
#   cat("Enter the ",i,"th group initial design \n")
#   cat("Enter the quantity of this design \n")
#   ni[i] <- as.integer(readline(prompt = ""))
#   ti=c()
#   for(j in 1:ni[i]){
#     tij <- as.numeric(readline(prompt = ""))
#     ti <- c(ti,tij)
#   }
#   
#   t <- c(t,ti) 
# }
t<-c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)

#Fixed effects parameters values
beta<-c(1.6,8,0.13)

#(Diagonal Matrix of) variance for inter-subject random effects:
omega<-diag(c(0.7,0.02,0.06))


#Random effect model (1) = additive  (2) = exponential 
#------------------------------------------------------------------
Trand<-2;

if ( Trand == 1 ) {
#   beta <- beta + rnorm(length(t)*length(paramName), mean = 0, sd = omega)
#   bm <- matrix(beta,nrow=3,ncol=12)
  form11 <- form1
  for(i in 1:length(paramName)){
    form11 <- gsub(paramName[i], paste0("(",paramName[i],"+b)"), form11)
    
  }
  
} else {
#   b = rnorm(length(t)*length(paramName), mean = 0, sd = omega)
#   beta <- beta * t(exp(b))
#   bm <- matrix(beta,nrow=3,ncol=12)
  form11 <- form1
  for(i in 1:length(paramName)){
    form11 <- gsub(paramName[i], paste0("(",paramName[i],"*exp(b))"), form11)
  }
}


#calculate matrix E for n individuals
equatf <- parse(text = form11, n=-1)
f<-function(paramF){eval(equatf[[1]])}

ka = beta[1]
V = beta[2]
Cl = beta[3]

param <- c(beta,t)
#calculate the observations with personnal parameters
b <- 0
fixed<-f(param)

#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
sig.inter<-0.6
sig.slope<-0.07

var<-diag((sig.inter+sig.slope*fixed)^2)

#calculate variance Var
form2<-paste0("( sig.inter + sig.slope * ", form11,")^2")
#form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
Vmodel<- parse(text = form2)

epsilon<-rnorm(length(t), mean = 0, sd = diag(var))

obs<-fixed+epsilon


df<-deriv(equatf[[1]],PSI)
mdf<-attributes(eval(df))$gradient
#delete the last two columns (correspond to sig.inter and sig.slope) 
mdfi <- mdf[,-c(length(PSI)-1,length(PSI))]
#calculate variance Vi
Vi <- mdfi %*% omega %*% t(mdfi) + var

#-------------------------------------delete----------
#calculate the matrix dVi/dlambda
# the first three columns of dVi/dlambda
dvi <- mdfi %*% diag(1,length(paramName)) %*% t(mdfi) 
#combine the last two columns of dVi
dv<-deriv(Vmodel[[1]],PSI)
mdv<-attributes(eval(dv))$gradient
dvi<-cbind(dvi,mdv[,length(PSI)-1],mdv[,length(PSI)])
dvi <- diag(diag(dvi))
# for(i in 1:length(paramName)){
#   mdv[,i]=0
# }
# m<- list()
# for(i in (length(paramName)+1):length(PSI)){
#   m[[i-length(paramName)]] <- diag(mdv[,i])
# }
#----------------------------------delete---------------

M_A <- t(mdfi) %*% solve(Vi) %*% mdfi 
for(i in 1:length(PSI)){
  M_A <- cbind(M_A,0)
  M_A <- rbind(M_A,0)
}

M_B <- matrix(rep(0),nrow=length(PSI)+length(paramName),ncol=length(PSI)+length(paramName))
for(i in 1:length(paramName)){
  
  for(j in 1:length(paramName)){
    M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfi[,i] %*% t(mdfi[,i])) %*% solve(Vi) %*% (mdfi[,j] %*% t(mdfi[,j])) %*% solve(Vi))))
  }
  
}

for(k in 4:5){
  for(l in 4:5){
    M_B[3+k,3+l] <- 1/2 * sum(diag(diag(mdv[,k]) %*% solve(Vi) %*% diag(mdv[,l]) %*% solve(Vi)))
  }
}
M_F <- M_A+M_B



