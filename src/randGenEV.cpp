#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]

/*
*************
* generator *
*************
input arguments:
* value: a real number between 0 and 1
* type: the type of function wanted 
* family: parametric family for the generators 
* parameter: a real parameter 
* type: 1 = derivative, 2 = generator, 3 = inverse
* family: 1 = frechet, 2 = cuadrasauge, 3 = sinus, 4 = exponential
*/

double generatorCPP(double value, int type, double parameter, int family){
  
    double f;
    double df;
    double invf;
    double out;
    
    if(family==1){
        f = (1-parameter)*value+parameter;
        df = 1-parameter;
        invf = (value-parameter)/(1-parameter);
      }else if(family==2){
        f = pow(value,1-parameter);
        df = (1-parameter)*pow(value,-parameter);
        invf = pow(value,1/(1-parameter));
      }else if(family==3){
        f = sin(value*parameter)/sin(parameter);
        df = parameter*cos(parameter*value)/sin(parameter);
        invf = asin(value*sin(parameter))/parameter;
      }else{//family==4
        f = exp((pow(value,parameter)-1)/parameter);
        df = exp((pow(value,parameter)-1)/parameter)*pow(value,parameter-1);
        invf = pow(parameter*log(value)+1,1/parameter);
      }
    
    if(type==1){
      out = df;
      }else if(type==2){
        out = f;
      }else{//type==3
        out = invf;
      }
    return out;
}

/*
*********************
* random generation *
*********************
*/
// [[Rcpp::export]]
NumericMatrix randGenCPP(int sampleSize, int family, NumericVector parameter,
NumericVector latent, NumericMatrix v){
  
  int n = sampleSize;
  int d = parameter.size();
  NumericMatrix data(n,d);
  double invf;
  double f;
  double df;
  
  for(int m=0;m<n;m++){
    for(int j=0;j<d;j++){
      invf = generatorCPP(v(m,j),3,parameter[j],family);
      df = generatorCPP(latent[m],1,parameter[j],family);
      f = generatorCPP(latent[m],2,parameter[j],family);
      if(v(m,j)>=f){
        data(m,j)=invf;
      }else if(v(m,j)<latent[m]*df){
        data(m,j)=v(m,j)/df;
      }else{
        data(m,j)=latent[m];
      }
    }
  }
  return data;
}

/*
***********************
* random generation 2 *
***********************
*/

// [[Rcpp::export]]
NumericMatrix randGenCPP_2(int sampleSize, int family, NumericVector parameter){
  
  int n = sampleSize;
  int d = parameter.size();
  NumericMatrix data(n,d);
  double invf;
  double f;
  double df;
  RNGScope scope; // I've seen that in Rcpp Gallery
  NumericVector tempLatent;
  double latent;
  NumericVector v(d);
  double tempV;
  
  for(int m=0;m<n;m++){
    tempLatent = runif(1);
    latent = tempLatent[0];
    v = runif(d);
    for(int j=0;j<d;j++){
      tempV = v[j];
      invf = generatorCPP(tempV,3,parameter[j],family);
      df = generatorCPP(latent,1,parameter[j],family);
      f = generatorCPP(latent,2,parameter[j],family);
      if(tempV>=f){
        data(m,j)=invf;
      }else if(tempV<latent*df){
        data(m,j)=tempV/df;
      }else{
        data(m,j)=latent;
      }
    }
  }
  return data;
}


/*
******************************************************
* random generation of the (extreme value) attractor *
******************************************************
*/

// [[Rcpp::export]]
NumericMatrix randGen_ev_CPP(int sampleSize, int family, NumericVector parameter,
NumericVector latent, NumericMatrix v, int underlyingSampleSize){
  int n1 = underlyingSampleSize;
  int n2 = sampleSize;
  int d = parameter.size();
  NumericMatrix data1(n1,d); // attracted data
  NumericMatrix data2(n2,d); // attracting extreme-value data
  NumericVector thisLatent(n1);
  NumericMatrix thisV(n1,d);
  double temp;
  
  for(int m=0;m<n2;m++){
    
    // provide the uniform numbers for randGen (thisLatent and thisV)
    for(int l=0;l<n1;l++){
      thisLatent[l] = latent[m*n1+l];
    }
    for(int b=0;b<d;b++){// runs the columns of thisV
      for(int a=0;a<n1;a++){// runs the rows of thisV
        thisV(a,b) = v(m*n1+a,b);
      }
    }
    
    // generate observations from the underlying copula
    data1 = randGenCPP(n1,family,parameter,thisLatent,thisV);
    
    for(int j=0; j<d;j++){// runs over the variables
      //compute the max of the column j
      temp=pow(data1(0,j),n1);
      for(int i=1;i<n1;i++){
        if(pow(data1(i,j),n1)>temp){
          temp=pow(data1(i,j),n1);
        }
      }
      data2(m,j)=temp; // the max of the column j is returned
    }
  }
  return data2;
}

/*
********************************************************
* random generation of the (extreme value) attractor 2 *
********************************************************
*/

// [[Rcpp::export]]
NumericMatrix randGen_ev_CPP_2(int sampleSize, int family, NumericVector parameter,
int underlyingSampleSize){
  int n1 = underlyingSampleSize;
  int n2 = sampleSize;
  int d = parameter.size();
  NumericMatrix data1(n1,d); // attracted data
  NumericMatrix data2(n2,d); // attracting extreme-value data
  double temp;
  
  for(int m=0;m<n2;m++){
        
    // generate observations from the underlying copula
    data1 = randGenCPP_2(n1,family,parameter);
    
    for(int j=0; j<d;j++){// runs over the variables
      //compute the max of the column j
      temp=pow(data1(0,j),n1);
      for(int i=1;i<n1;i++){
        if(pow(data1(i,j),n1)>temp){
          temp=pow(data1(i,j),n1);
        }
      }
      data2(m,j)=temp; // the max of the column j is returned
    }
  }
  return data2;
}
