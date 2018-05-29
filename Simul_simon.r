# Simulation with sigmoid Emax model

###################################################################
# PD function : sigmoid Emax model
###################################################################

rm(list=ls(all=TRUE))
model<-function(VA,k,be,Emax,ED50,dose){
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
  bi<-rmvnorm(Nsub, mean = c(0,0,0,0,0), sigma=BSV)

  mu<-matrix(rep(c(VA,k,be,Emax,ED50),Nsub),byrow=T,ncol=4)

  pat<-rep(1:Nsub)

  thetaij<-mu*exp(bi)
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

simul_effectE<-function(Nsub,dose,a_sig,seed,Directory1,Directory2){
  setwd(Directory1)
  thetaij<-read.table(paste("param",seed,".txt",sep=""),header=T,sep="\t")
  E<-c()
  for(i in 1:(Nsub)){
  E<-c(E,model(thetaij$VA[i],thetaij$k[i],thetaij$be[i],thetaij$Emax[i],thetaij$ED50[i],dose))
  }
  set.seed(seed)
  Error<-rnorm(Nsub*length(dose),mean=0,sd=a_sig)
  E_sim<-E + Error
  E<-E_sim
  pat<-rep(1:Nsub,each=length(dose))

  Dose<-rep(dose,times=Nsub)
  Data<-as.data.frame(cbind(pat,Dose,E_sim,E))
  names(Data)<-c("subject","dose","E_sim","E")
  setwd(Directory2)
  write.table(Data,paste("data",seed,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
  rm(Data)
  }

simul_effectE2<-function(Nsim,Nsub,dose,sig,seed,Directory1,Directory2){
  for(i in 1:Nsim){
  simul_effectE(Nsub,dose,sig,i,Directory1,Directory2)
  }
}

###################################################################
# Simulation of dose-response trial data for MONOLIX
###################################################################
format_monolix<-function(Nsub,NsampM,seed,Directory1,Directory2){
  setwd(Directory1)
  Data<-read.table(paste("data",seed,".txt",sep=""),header=T,sep="\t")
  Data<-Data[,-3]
  rownames(Data)<-NULL
  Data<-as.data.frame(Data)
  names(Data)<-c("ID","X","Y")

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
k_s<-0.0005
be_s<-0.2
Emax_s<-30
ED50_s<-150


# variability
BSV_s<-diag(c(0.07,0.5,1,150,0))
a_sig_s<-sqrt(28)

# doses
dose_s<-c(0, 25, 50, 100, 150, 300, 500)

# number of subjects, of samples, of simulations
Nsub_s<-300
Nsamp_s<-length(dose_s)
NsampM_s<-Nsamp_s
Nsim_s<-1000

# graine
seed_s<-1


# Simulation for gamma = 3, 4 doses, 1000 dataset of 100 subjects
###################################################################
g_s<-3
dose_s<-c(0,100,300,1000)
Nsub_s<-100
Nsamp_s<-length(dose_s)
NsampM_s<-Nsamp_s
Nsim_s<-1000

# Directories
Directory1_s<-"C:\\Users\\Thu Thuy Nguyen\\Desktop\\Thu Thuy\\TravauxThèse\\CalculMF\\Simulation\\Design4doses_gamma3\\Param"
Directory2_s<-"C:\\Users\\Thu Thuy Nguyen\\Desktop\\Thu Thuy\\TravauxThèse\\CalculMF\\Simulation\\Design4doses_gamma3\\Data_original"
Directory3_s<-"C:\\Users\\Thu Thuy Nguyen\\Desktop\\Thu Thuy\\TravauxThèse\\CalculMF\\Simulation\\Design4doses_gamma3\\Data_Monolix"

# appel des fonctions
simul_param(Nsub_s,Emax_s,ED50_s,g_s,E0_s,BSV_s,seed_s,Directory1_s)
simul_param2(Nsim_s,Nsub_s,Emax_s,ED50_s,g_s,E0_s,BSV_s,seed_s,Directory1_s)

simul_effectE(Nsub_s,dose_s,a_sig_s,seed_s,Directory1_s,Directory2_s)
simul_effectE2(Nsim_s,Nsub_s,dose_s,a_sig_s,seed_s,Directory1_s,Directory2_s)

format_monolix(Nsub_s,NsampM_s,seed_s,Directory2_s,Directory3_s)
format_monolix2(Nsim_s,Nsub_s,NsampM_s,Directory2_s,Directory3_s)





