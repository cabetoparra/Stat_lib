#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

void Loren_dens(double *x,double *y,double gamma,double *eval,int dim,int size);
double mtwert(double *data,int length);
void stupidsort(double *a, int length);
void Loren_unfold(double gamma,double *eval,int dim,int size,double *mls,double *err,double *emax);
void simp(int n,double h,double f[],double *s);
void hist(double *array,int size,double xi,double xf,int bin,char *file);
void easyunfold(double gamma,double *eval,int dim,double *mls,double *err,double *emax);
