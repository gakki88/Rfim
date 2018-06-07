getEstimation<-function(N){
  
  directory <-paste0("D:\\yuxin\\Simulation\\resultatCDC\\md_",1,"\\projectSimulSimon")
  setwd(directory)
  Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
  datas<-Data[,2]
  
  for(i in 2:N){
    tryCatch({
      directory <-paste0("D:\\yuxin\\Simulation\\resultatCDC\\md_",i,"\\projectSimulSimon")
      setwd(directory)
      Data<-read.table(paste("estimates.txt",sep=""),header=T,sep=";")
      datas<-cbind(datas,Data[,2])
    },error=function(e){
      
    }
      )
    
  }
    return(datas)
}

gg<-getEstimation(500)
SEemp<-c()
for(i in 1:10){
  SEemp <- c(SEemp,sqrt(var(gg[i,])))
  boxplot(gg[i,])
}

SEemp[1:5]/c(55,0.005,0.2,30,150)*100
hist(gg[1,])
dim(gg)
