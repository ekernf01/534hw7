
#include "matrices.h"
#include <stdio.h>
#include <iostream>
void eric_gsl_mat_print(char* outputfilename, gsl_matrix* mat_to_print);
double** readmatrix(char* filename,int n,int p);
void freematrix(int n,double** m);
gsl_matrix* eric_covariance(double** data, int num_vars, int num_cases);
