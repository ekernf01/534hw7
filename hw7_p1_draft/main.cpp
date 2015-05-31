/*
 FILE: MAIN.CPP

 This program draws 10,000 multivariate Normal variates with (theoretical)
 covariance matching the sample covariance from erdata.txt.
 
*/


#include "matrices.h"
#include "regmodels.h"
#include <stdio.h>
#include <iostream>
using namespace std;

int main()
{
  int n = 158; //sample size
  int p = 51; //number of variables
  int num_preds = p-1;//number of predictors
  char datafilename[] = "erdata.txt"; //name of the data file
  char outputfilename[] = "eric_cov_10_000.txt";

  //allocate the data matrix
  double** data = allocmatrix(n,p);

  //read the data
  readmatrix(datafilename,n,p,data);
  
  //Discard the first column, which contains the response
  int* row_inds = new int[n];
  int* col_inds = new int[num_preds];
  for(int i=1; i++; i<=n){
    row_inds[i] = i;
  }
  for(int i=1; i++; i<=(num_preds)){
    col_inds[i] = i;
  }
  gsl_matrix* predictors_only = MakeSubmatrix(data,n, row_inds, num_preds, col_inds);

  
  
  //Obtain the covariance matrix
  gsl_matrix* cov_or_sqrt_cov = gsl_matrix_alloc(num_preds,num_preds)
  for(int i=1; i++; i<=num_preds){
    for(int j=1; j++; j<=num_preds){
      if(i<=j){
        cov_entry = gsl_stats_covariance(col_array[i], 1, col_array[j], 1, n); //The extra ones say not to skip elements.
        gsl_matrix_set(cov_or_sqrt_cov, i-1, j-1, cov_entry);//GSL indexes from 0
        gsl_matrix_set(cov_or_sqrt_cov, j-1, i-1, cov_entry);
      }
    }
  }

  //free memory
  freematrix(n,data);
  gsl_matrix_free(predictors_only);
  delete col_inds;
  delete row_inds;
  
  return(1);
}
