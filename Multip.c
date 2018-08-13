#include <stdio.h>
#include "matrix.h"




float mutilicative(float *MFi, float w,int iteration, float lambda, float delta)
{
  float d;
  int it=0;
  for(int k=1; k < iteration ; k++){
  
  /*Calculate the sum of matrix with weight associated
  idea:   Mw<-Mw+w[i]*MFi[[i]]*/
    MFw<-Map("*",MFi,w)     #associate weights
      Mw<-Reduce("+",MFw)     #sum of all matrices
      dm<-dim(Mw)[1]
    Dphi <- det(Mw)^(1/dm) * solve(Mw)/dm            # calculate derivatives of function phi_D
      d<-sapply(MFi, FUN = function(x) sum(diag(Dphi %*% x)))      #the vector of multiplier
      w<- w*d^lambda/sum(w*d^lambda)          #develop weight
      
      if(max(d)<(1+delta)*sum(w*d)){break}    #stop criterion
        it<-it+1                                #number of iteration calculator
      if(it%%100==0) print(it)
  }
  return(w);
}