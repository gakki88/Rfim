funFIMsm <- function(equation,parameter,Vparam,t_indiv,nbSubjects){
#Name of the fixed effects parameters
paramName<-c("ka","V","Cl")

paramF<-c(paramName,"t")

#model equation
#form1<- paste0("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))")
#equatf <- parse(text = form1, n=-1)
equatf <- parse(text = equation, n=-1)
f<-function(paramF){eval(equatf[[1]])}

form2<-paste0("( sig.inter + sig.slope * ", form1,")^2")
#form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
Vmodel<- parse(text = form2)

PSI<-c(paramName,"sig.inter", "sig.slope")


#--------------------------------------------
#calculate Fisher Information Matrix
#--------------------------------------------
dose<-c(100)

#Fixed effects parameters values
#beta<-c(1.6,8,0.13)
beta <- parameter
ka=beta[1]
V=beta[2]
Cl=beta[3]

#sample times for an individual design
#t<-c(0.5, 2, 9, 24, 48, 96)
t <- t_indiv

#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
#sig <- c(0.6,0.07)
sig <- Vparam
sig.inter<-sig[1]
sig.slope<-sig[2]


#lambda<-c(sig.inter, sig.slope)

#psi<-c(beta,lambda)
param <- c(beta,t)
fixed<-f(param)
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

deterFim <- det(M_F)
SE <- diag(sqrt(solve(M_F)))
RSE <- SE / diag(M_F) *100

return(M_F*nbSubjects)
}

funFIMsm("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c(1.6,8,0.13),c(0.6,0.07),c(0.5, 2, 9, 24, 48, 96),16)





#Diagonal Matrix of variance for inter-subject random effects:
# omega<-diag(c(0.7,0.02,0.06))
# 
# 
# 
# ## Higher derivatives:
# DD <- function(expr, name, order = 1) {
#   if(order < 1) stop("'order' must be >= 1")
#   if(order == 1) D(expr, name)
#   else DD(D(expr, name), name, order - 1)
# }
# DD(expression(sin(x^2)), "x", 3)
# 
# func2 <- function(x) c(sin(x), cos(x))
# x <- (0:1)*2*pi
# jacobian(func2, x)
# jacobian(func2, x, "complex")
