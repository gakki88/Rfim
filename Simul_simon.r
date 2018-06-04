# Simulation with sigmoid Emax model

###################################################################
# PD function : sigmoid Emax model
###################################################################

rm(list=ls(all=TRUE))
model<-function(VA,k,be,Emax,ED50,dose,t){
  f<-VA+(1-exp(-k*t))*(Emax*dose/(ED50+dose)-be*VA)
  return(f)
}


######################################################################
#  Simulation of parameters for sigmoid Emax model
######################################################################

simul_param<-function(Nsub,VA,k,be,Emax,ED50,BSV,seed,Directory){
  library(nlme)
  library(mvtnorm)
  
  set.seed(seed)
  bi<-rmvnorm(sum(Nsub), mean = c(0,0,0,0,0), sigma=BSV)
  
  mu<-matrix(rep(c(VA,k,be,Emax,ED50),sum(Nsub)),byrow=T,ncol=5)
  
  pat<-rep(1:sum(Nsub))
  
  thetaij<-mu*exp(bi)
  thetaij[,4]<-mu[,4]+bi[,4]
  thetaij<-as.data.frame(cbind(pat,thetaij))
  names(thetaij)<-c("subject","VA","k","be","Emax","ED50")
  setwd(Directory)
  write.table(thetaij,paste("param",seed,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
  rm(thetaij)
  
}

simul_param2<-function(Nsim,Nsub,VA,k,be,Emax,ED50,BSV,seed,Directory){
  for(i in 1:Nsim){
    simul_param(Nsub,VA,k,be,Emax,ED50,BSV,i,Directory)
  }
}


#########################################################################
# Simulations of responses for sigmoid Emax model
#########################################################################

simul_effectE<-function(Nsub,dose,t,a_sig,seed,Directory1,Directory2){
  setwd(Directory1)
  thetaij<-read.table(paste("param",seed,".txt",sep=""),header=T,sep="\t")
  E<-c()
  for(j in 1:length(Nsub)){
    for(i in 1:(Nsub[j])){
      E<-c(E,model(thetaij$VA[i],thetaij$k[i],thetaij$be[i],thetaij$Emax[i],thetaij$ED50[i],dose[j],t[[j]]))
    }
  }
  
  set.seed(seed)
  # 4 for the number of sampling times correspond to each dose
  Error<-rnorm(Nsub*length(dose)*4,mean=0,sd=a_sig)
  E_sim<-E + Error
  E<-E_sim
  pat<-rep(1:sum(Nsub),each=4)
  
  Dose<-c()
  Time<-c()
  for(j in 1:length(t)){
    Dose<-c(Dose,rep(dose,times=length(Nsub),each=Nsub[j]*4))
    Time<-c(Time,rep(t[[j]],times=Nsub[j]))
  }
  Data<-as.data.frame(cbind(pat,Dose,Time,E_sim))
  names(Data)<-c("ID","Dose","Time","Y")
  setwd(Directory2)
  write.table(Data,paste("data",seed,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
  rm(Data)
}

simul_effectE2<-function(Nsim,Nsub,dose,t,sig,seed,Directory1,Directory2){
  for(i in 1:Nsim){
    simul_effectE(Nsub,dose,t,sig,i,Directory1,Directory2)
  }
}

###################################################################
# Simulation of dose-response trial data for MONOLIX
###################################################################
format_monolix<-function(Nsub,NsampM,seed,Directory1,Directory2){
  setwd(Directory1)
  Data<-read.table(paste("data",seed,".txt",sep=""),header=T,sep="\t")
  Data<-Data[,-4]
  rownames(Data)<-NULL
  Data<-as.data.frame(Data)
  names(Data)<-c("ID","Dose","Time","Y")
  
  setwd(Directory2)
  write.table(Data,paste("DATA",seed,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
}


format_monolix2<-function(Nsim,Nsub,NsampM,Directory1,Directory2){
  for(i in 1:Nsim){
    format_monolix(Nsub,NsampM,i,Directory1,Directory2)
  }
}


###################################################################
# Simulation sigmoid Emax model Execution
###################################################################

# PD parameters
VA_s<-55
k_s<-0.005
be_s<-0.2
Emax_s<-30
ED50_s<-150


# variability
BSV_s<-diag(c(0.07,0.5,1,150,0))
a_sig_s<-sqrt(28)

# doses
dose_s<-c(0, 50, 100, 500)

# sampling times
t_s<-list(c(0, 168, 196, 672),c(0, 7, 644, 672),c(0, 7, 644, 672),c(0, 168, 196, 672))

# number of subjects, of samples, of simulations
#Nsub_s<-c(75,75,75,75)
Nsub_s<-c(107,47,14,132)
Nsamp_s<-length(dose_s)
NsampM_s<-Nsamp_s
Nsim_s<-500

# graine
seed_s<-1


# Simulation for gamma = 3, 4 doses, 1000 dataset of 100 subjects
###################################################################
# g_s<-3
# dose_s<-c(0,25,50,100,150,300,500)
# Nsub_s<-c(75,75,75,75)
# Nsamp_s<-length(dose_s)
# NsampM_s<-Nsamp_s
# Nsim_s<-1000

# Directories
# Directory1_s<-"D:\\yuxin\\Simulation\\Param"
# Directory2_s<-"D:\\yuxin\\Simulation\\Data_original"
Directory1_s<-"D:\\insa\\inserm\\Simulation\\Param"
Directory2_s<-"D:\\insa\\inserm\\Simulation\\Data_original"
#Directory3_s<-"D:\\yuxin\\Simulation\\Data_Monolix"

# appel des fonctions
simul_param(Nsub_s,VA_s,k_s,be_s,Emax_s,ED50_s,BSV_s,seed_s,Directory1_s)
#simul_param2(Nsim_s,Nsub_s,VA_s,k_s,be_s,Emax_s,ED50_s,BSV_s,seed_s,Directory1_s)

simul_effectE(Nsub_s,dose_s,t_s,a_sig_s,seed_s,Directory1_s,Directory2_s)
#simul_effectE2(Nsim_s,Nsub_s,dose_s,t_s,a_sig_s,seed_s,Directory1_s,Directory2_s)

#format_monolix(Nsub_s,NsampM_s,seed_s,Directory2_s,Directory3_s)
#format_monolix2(Nsim_s,Nsub_s,NsampM_s,Directory2_s,Directory3_s)

# Plot spaghetti graphique
##################################################################
setwd("D:\\insa\\inserm\\Simulation\\Data_original")
Data<-read.table(paste("data",1,".txt",sep=""),header=T,sep="\t")
Data1<-Data[which(Data$Dose==500),]
plot(x=Data1$Time,y=Data1$Y,col="grey50")
lines(x=Data1$Time,y=Data1$Y,col="grey80")

# Plot the prediction line 
t<-seq(0,672,1)
equatf <- parse(text = "VA_s+(1-exp(-k_s*t))*(Emax_s*dose/(ED50_s+dose)-be_s*VA_s)", n=-1)
f<-function(paramF){eval(equatf[[1]])}
dose<-500
x<-c(VA_s,k_s,be_s,Emax_s,ED50_s,dose,t)
fvalue<-f(x)
# the line of the model function with manually ajusted y_axis
lines(t,fvalue,type="l",lwd=3) 
