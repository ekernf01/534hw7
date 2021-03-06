#include "hw7p1.h"
#include <stdio.h>
using namespace std;

//Prints a GSL matrix to the file outputfilename
void eric_gsl_mat_print(char* outputfilename, gsl_matrix* mat_to_print){
  char format[] = "%f";
  FILE* out = fopen(outputfilename,"w");
  if(NULL==out){
    printf("Cannot open output file [%s]\n",outputfilename);
    exit(1);
  }
  gsl_matrix_fprintf(out, mat_to_print, format);
  fclose(out);
  return;
}


//reads from a file a matrix with n rows and p columns
//allocates and returns the 2d array m so that m[p-1][n-1] is valid
double** readmatrix(char* filename,int n,int p){
  double** m = new double*[p];
  int i,j;
  double s;
  FILE* in = fopen(filename,"r");
  
  if(NULL==in){
    printf("Cannot open input file [%s]\n",filename);
    exit(1);
  }
  for(j=0;j<p;j++){
    m[j] = new double[n];
    for(i=0;i<n;i++){
      fscanf(in,"%lf",&s);
      m[j][i] = s;
    }
  }
  fclose(in);
  return m;
}

//frees the memory for a matrix with n rows
void freematrix(int n,double** m){
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
	return;
}

//This function returns a GSL matrix containing the covariance of the variables stored in data.
//data[num_vars-1][num_cases-1] should return a double.
//This function allocates memory.
gsl_matrix* eric_covariance(double** data, int num_vars, int num_cases){
  gsl_matrix* cov = gsl_matrix_alloc(num_vars,num_vars);
  double cov_entry = 0;
  for(int i=0; i<num_vars; i++){
    for(int j=0; j<num_vars; j++){
      if(i<=j){
        cov_entry = gsl_stats_covariance(data[i], 1, data[j], 1, num_cases); //The extra ones say not to skip elements.
        gsl_matrix_set(cov, i, j, cov_entry);
        gsl_matrix_set(cov, j, i, cov_entry);
      }
    }
  }
  return cov;
}