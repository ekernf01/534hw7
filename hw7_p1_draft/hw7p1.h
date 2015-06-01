
#include "matrices.h"
#include <stdio.h>
#include <iostream>
double** readmatrix(char* filename,int n,int p);
void freematrix(int n,double** m);
gsl_matrix* eric_covariance(double** data, int num_vars, int num_cases);
