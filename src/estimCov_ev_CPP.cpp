
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix estimCov_ev_CPP(NumericMatrix data){
  int n = data.nrow();
  int d = data.ncol();
  int p = d*(d-1)/2;
  int pairIndex1;
  int pairIndex2;
  NumericMatrix covMatLambdaHat(p,p);
  double pseudoLambdaij;
  double pseudoLambdakl;
  double lambda_ij;
  double lambda_kl;
  double etaijkl;
  double temp;
  NumericVector theMaxij(n);
  NumericVector theMaxkl(n);

  pairIndex1 = 0;
  for(int i=0;i<d-1;i++){
    for(int j=i+1;j<d;j++){
      pairIndex1 = pairIndex1 + 1;
      
      
      // computation of lambda_ij
      for(int m=0;m<n;m++){
        if(data(m,i)>=data(m,j)){
          theMaxij[m]=data(m,i);
        }else{
          theMaxij[m]=data(m,j);
        }
      }
      temp = 0;
      for(int m=0;m<n;m++){
        temp = temp + theMaxij[m];
      }
      pseudoLambdaij = temp/n;
      lambda_ij = 3-1/(1-pseudoLambdaij);
      
      pairIndex2 = 0;
      for(int k=0;k<d-1;k++){
        for(int l=k+1;l<d;l++){
          pairIndex2 = pairIndex2 + 1;
          
          // computation of lambda_kl
          for(int m=0;m<n;m++){
            if(data(m,k)>=data(m,l)){
              theMaxkl[m]=data(m,k);
            }else{
              theMaxkl[m]=data(m,l);
            }
          }      
          temp = 0;
          for(int m=0;m<n;m++){
            temp = temp + theMaxkl[m];
          }
          pseudoLambdakl = temp/n;
          lambda_kl = 3-1/(1-pseudoLambdakl);
          
          // computation of etaijkl
          temp = 0;
          for(int m=0;m<n;m++){
            temp = temp + theMaxkl[m]*theMaxij[m];
          }
          etaijkl=temp/n;
          //printf("\n%d et %d",pairIndex1,pairIndex2);
          covMatLambdaHat(pairIndex1-1,pairIndex2-1)=pow(3-lambda_ij,2)*pow(3-lambda_kl,2)* ( etaijkl - (  (2-lambda_ij)*(2-lambda_kl) ) / ( (3-lambda_ij)*(3-lambda_kl) ) );    
        }
      }
    }
  }
  return covMatLambdaHat;
}
