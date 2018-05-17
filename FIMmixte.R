funFIMem <- function(equation,paramName,beta,o,sigma,t_group,Trand,d,PropSubjects,nbTot){
  #Name of the fixed effects parameters and sampling times
  paramF<-c(paramName,"t")
  
  #model equation
  form1<- equation
  
  PSI<-c(paramName,"sig.inter","sig.slope")
  lpsi<-length(PSI)

  
  #(Diagonal Matrix of) variance for inter-subject random effects:
  omega<-diag(o)
  
  

  
  #gather all groups of protocol
  M_f<-list()
  M_F <- matrix(rep(0),nrow=lpsi+length(paramName),ncol=lpsi+length(paramName))
  for(q in 1:length(t_group)){
    t<-c()
    t<-c(t_group[[q]])
    
    #dose value
    dose<-d[q]
    
    #calculate matrix E for n individuals
    equatf <- parse(text = form1, n=-1)
    f<-function(paramF){eval(equatf[[1]])}
    
    #Fixed effects parameters values
    for(i in 1:length(paramName)){
      assign(paramName[i],beta[i])
    }
    
    #calculate the observations with personnal parameters
    param <- c(beta,t)
    fixed<-f(param)
    
    #Standard deviation of residual error (sig.inter+sig.slope*f)^2:
    sig.inter<-sigma[1]
    sig.slope<-sigma[2]
    
    var<-diag((sig.inter+sig.slope*fixed)^2)
    
    
    #calculate variance Var
    form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
    Vmodel<- parse(text = form2)
    
    #get derivatives for fixed parameters
    df<-deriv(equatf[[1]],PSI)
    mdf<-attributes(eval(df))$gradient
    #delete the last two columns (correspond to sig.inter and sig.slope) 
    mdfi <- mdf[,-c(length(PSI)-1,length(PSI))]
    #complete derivative for exponential random effect model   
    #Random effect model (1) = additive  (2) = exponential 
    #------------------------------------------------------------------
    beta0 <- beta
    beta0[which(Trand==1)] <- 1 
    #if(Trand ==2 ){
      mdfie <- mdfi %*% diag(beta0)
    #}else {mdfie <- mdfi}
     
    #calculate variance Vi
    Vi <- mdfie %*% omega %*% t(mdfie) + var
    
    #get derivatives of sigma
    dv<-deriv(Vmodel[[1]],PSI)
    mdv<-attributes(eval(dv))$gradient
    
    
    #calculate matrix part A
    M_A <- t(mdfi) %*% solve(Vi) %*% mdfi
     
    #complete the rest of the matrix with 0
    for(i in 1:length(PSI)){
      M_A <- cbind(M_A,0)
      M_A <- rbind(M_A,0)
    }
    
    #calculate matrix part B
    #initialize the matrix with 0
    M_B <- matrix(rep(0),nrow=length(PSI)+length(paramName),ncol=length(PSI)+length(paramName))
    #prepare a list of matrix of derivatives of sigma to simplify usage
    m<-list()
    for(i in (length(paramName)+1):lpsi){
      if(length(t)==1){
        m[[i]] <- mdv[i]
      }else{
        m[[i]] <- diag(mdv[,i])
      }
      
    }
    #calculate first three rows of part B
    for(i in 1:length(paramName)){
      
      for(j in 1:length(paramName)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi))))
      }
      for(j in (length(PSI)-1):length(PSI)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% m[[j]] %*% solve(Vi))))
      }
      
    }
    #calculate the last two rows of partB
    for(i in (length(PSI)-1):length(PSI)){
      for(j in 1:length(paramName)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag( m[[i]] %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi)))
      }
      for(j in (length(PSI)-1):length(PSI)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(m[[i]] %*% solve(Vi) %*% m[[j]] %*% solve(Vi)))
      }
    }
    
    M_f[[q]] <- (M_A+M_B)*PropSubjects[q]
    M_F <-M_F+M_f[[q]]
  }
  M_F <- M_F *nbTot
  
  #set names for vectors 
  fname<-c()
  for(n in 1:length(paramName)){
    fname<-c(fname,paste0("u_",paramName[n]))
  }
  for(n in 1:length(paramName)){
    fname<-c(fname,paste0("w2_",paramName[n]))
  }
  fname<-c(fname,"sig.inter","sig.slope")
  rownames(M_F) <- fname
  colnames(M_F) <- fname
  
  if(0 %in% c(o,sigma)){
    M_F<-M_F[-c(length(paramName)+which(c(o,sigma)==0)),-c(length(paramName)+which(c(o,sigma)==0))]
    PSI<- PSI[-c(length(paramName)+which(c(o,sigma)==0))]
  }

  if(length(t)==1 || (cor(M_F)==1) ){
    return(list(M_F,det(M_F)))
  }else{
    deterFim <- det(M_F)
    SE <- sqrt(diag(solve(M_F)))
    RSE <- 100 * SE / c(beta,o,sigma)[which(c(beta,o,sigma)!=0)]
    CritereDopt <- deterFim^(1/(length(PSI)+length(paramName)))
    

#     t<-seq(min(unlist(t_group)),max(unlist(t_group)),(max(unlist(t_group))-min(unlist(t_group)))/100)
#     x<-c(beta,t)
#     fvalue<-f(x)
#     lines(t,fvalue) #then plot the line of function
    return(list(M_F,deterFim,CritereDopt,SE,RSE))
    #return(M_F)
    #write the output in console
    
    # cat("******************* FISHER INFORMATION MATRIX ******************\n")
    # print(M_F)
    # 
    # cat("\n\n******************* DETERMINANT OF THE MATRIX ******************\n", deterFim,
    #     "\n\n******************* STANDARD ERROR ******************\n",SE,
    #     "\n\n******************* RELATIVE STANDARD ERROR ******************\n",RSE,
    #     "\n\n******************* CRITERION ******************\n",CritereDopt)
    
  }
  
  
  
  
}
#exercice1
#funFIMem("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.7,0.02,0.06),c(0.6,0.07),list(c(0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120)),2,c(100),c(1),32)

#exercice2
#funFIMem("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))",c("ka","V","Cl"),c(1.6,8,0.13),c(0.7,0.02,0.06),c(0.6,0.07),list(c(0.5, 2, 9, 24, 48, 96),c(1, 6, 12, 36, 72, 120)),2,c(100,100),c(0.5,0.5),32)

#exemple dose-reponse
#funFIMem("Emax*t/(t+C50)+S0",c("Emax","C50","S0"),c(30,500,5),c(0.09,0.09,0.09),c(1,0),list(c(0,100 , 300 , 500 , 1000 , 2500 , 5000)),2,c(0),c(1),1)
