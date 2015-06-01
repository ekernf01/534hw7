/*
 FILE: MAIN.CPP

 This program draws 10,000 multivariate Normal variates with (theoretical)
 covariance matching the sample covariance from erdata.txt.
 It is intended to be readable rather than maximally efficient.
 One stumbling block: GSL's column and row-getting functions return GSL vectors,
 but GSL's covariance function demands simple arrays, as in double*.
 That mismatch prompted much of the bad design that went into this program.
 
 I have taken some of Adrian's code for (non-GSL) matrices.
 But, I did not put any of it in matrices.cpp. It's in hw7p1.cpp.
 
*/


#include "matrices.h"
#include <stdio.h>
#include <iostream>
#include "hw7p1.h"

using namespace std;

int main()
{
  int n = 158; //sample size
  int p = 51; //number of variables
  char datafilename[] = "erdata.txt"; //name of the data file
  char outputfilename[] = "eric_cov_10_000.txt";
  
  cout << "done reading data" << endl;
  
  //read the data
  double** data = readmatrix(datafilename,n,p);

  //Get the covariance
  gsl_matrix* cov_or_sqrt_cov = eric_covariance(data, p, n);
  
  //Do Cholesky. GSL leaves garbage above the diagonal. The loops below remove it.
  gsl_linalg_cholesky_decomp(cov_or_sqrt_cov);
  for(int i=0; i<p; i++){
    for(int j=0; j<p; j++){
      if(i<j){
        gsl_matrix_set(cov_or_sqrt_cov, i, j, 0);
      }
    }
  }

  //rng setup
  const gsl_rng_type * my_gsl_rng_type;
  gsl_rng * my_gsl_rng;
  gsl_rng_env_setup();
  my_gsl_rng_type = gsl_rng_default;
  my_gsl_rng = gsl_rng_alloc(my_gsl_rng_type);
  
  //Generate the actual random numbers
  int samp_size = 10000;
  gsl_matrix* sample_isotropic = gsl_matrix_alloc(samp_size, p);
  for(int i=0; i<samp_size; i++){
    for(int j=0; j<p; j++){
      gsl_matrix_set(sample_isotropic, i, j, gsl_ran_ugaussian(my_gsl_rng));
    }
  }
  
  //Transform them to have the right covariance
  //copy them, transposed, into a crude array
  gsl_matrix* sample_correct = gsl_matrix_alloc(samp_size, p);
  matrixproduct(cov_or_sqrt_cov, sample_isotropic, sample_correct);
  double** sample_correct_array = new double*[p];
  for(int i=0; i<p; i++){
    sample_correct_array[i] = new double[samp_size];
    for(int j=0; j<samp_size; j++){
      sample_correct_array[i][j] = gsl_matrix_get(sample_correct, j,i);
    }
  }
  
  //get the sample covariance for the synthetic data
  gsl_matrix* synth_samp_cov = eric_covariance(sample_correct_array, p, n);

  //print it
  char format[] = "%f";
  FILE* out = fopen(outputfilename,"w");
  if(NULL==out){
    printf("Cannot open output file [%s]\n",outputfilename);
    exit(1);
  }
  gsl_matrix_fprintf (out, synth_samp_cov, format);
  fclose(out);
  
  //free memory
  freematrix(p,data);
  gsl_matrix_free(synth_samp_cov);
  freematrix(p,sample_correct_array);
  gsl_matrix_free(sample_isotropic);
  gsl_matrix_free(sample_correct);
  
  return(1);
}
