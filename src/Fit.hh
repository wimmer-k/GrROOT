#ifndef FIT
#define FIT
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int Fit(double *x, double *y, double *sig, double *a, int MA, int NPoints, double *chisq); 
void gaussian(double x, double a[], double *y, double dyda[], int na);
void f1(double x[], double y[], double sig[], int ndata, double a[], int ia[],
       int ma, double **covar, double **alpha, double *chisq,
       void (*funcs)(double, double[], double*, double[], int),
       double *alamda);
void f2(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int));
void f3(double **a, int n, double **b, int m);
void f4(double **covar, int ma, int ia[], int mfit);
void f5(const char* error_text);
int *f6(long nl, long nh);
double *f7(long nl, long nh);
double **f8(long nrl, long nrh, long ncl, long nch);
double **f9(long nrl, long nrh, long ncl, long nch);
void ff6(int *v, long nl, long nh);
void ff7(double *v, long nl, long nh);
void ff8(double **m, long nrl, long nrh, long ncl, long nch);
void ff9(double **m, long nrl, long nrh, long ncl, long nch);

#endif
