#include <stdlib.h>
#include <iostream>
using namespace std;
#include <R.h>

/* Some very simple linear algebra functions */

/* If resid==0, return a vector containing the appropriate stratum (weighted) 
   means. Otherwise, center the input  vector around these. i.e  calculate
   either the "fitted value" or the residual from a model in which only 
   strata are fitted. In this and following functions, ynew can coincide with 
   y. Matrices are stored in Fortran order  
   Returns number of empty strata */

extern "C" {

  int wcenter(const double *y, int n, const double *weight, const int *stratum, int nstrata, int resid, double *ynew) {
    int i = 0, s=0;
    if (!stratum) {
      if (!nstrata) {
	/* Nothing to do ... if necessary copy input to output */
	if (ynew!=y) 
	  for(i=0; i<n; i++) ynew[i]  = resid? y[i]: 0.0; /* 0.0 ??? */
	return(0);
      }
      else
	nstrata = 1;
    }
    int empty = 0;
  if (nstrata>1) {
    double *swy, *swt;
    swy = (double *) Calloc(nstrata, double);
    swt =  (double *) Calloc(nstrata, double);
    for (s=0; s<nstrata; s++) 
      swy[s] = swt[s] = 0.0;
    if (weight) {
      for (i=0; i<n; i++) {
	double wi  = weight[i];
	int s = stratum[i]-1;
	swt[s] += wi;
	swy[s] += wi*y[i];
      }
    }
    else {
      for (i=0; i<n; i++) {
	int s = stratum[i]-1;
	swt[s] ++;
	swy[s] += y[i];
      }
    }
    for (s=0; s<nstrata; s++) {
      double sws = swt[s];
      if (sws > 0.0) 
	swy[s] /= sws;
      else
	empty++;
    }
    for (i=0; i<n; i++) {
      int s = stratum[i] -1;
      if (swt[s]) 
	ynew[i] = resid? y[i] - swy[s]: swy[s];
    }
    Free(swy);
    Free(swt);
  }
  else {    //only one stratum?
    double swt=0.0, swy=0.0;
    if (weight) {
      for (i=0; i<n; i++) {
	double wi = weight[i];
	swt += wi;
	swy += wi*y[i];
      }
    }
    else {
      for (i=0; i<n; i++) {
	swy += y[i];
      }
      swt = (double) n;
    }
    swy /= swt;
    if (swt>0) {
      for (i=0; i<n; i++){ynew[i] = resid? y[i] - swy: swy;}
    }
    else {empty = 1;}
  }
  return(empty);
}

/* Replace y by residual from (weighted) regression through the origin */

int wresid(const double *y, int n, const double *weight, const double *x, 
	   double *ynew) {
  double  swxx, swxy;
  swxy = swxx = 0.0;
  int i;
  if (weight) {
    for (i=0; i<n; i++) {
      double wi = weight[i];
      double xi = x[i];
      double wx = wi*xi;
      swxy += wx*y[i];
      swxx += wx*xi;
    }
  }
  else {
    for (i=0; i<n; i++) {
      double xi = x[i];
      swxy += xi*y[i];
      swxx += xi*xi;
    }
  }
  if (swxx>0) {
    swxy /= swxx;
    for (i=0; i<n; i++) {
      if (weight[i] > 0.) {ynew[i] = y[i] - swxy*x[i];} else {ynew[i] = y[i];}
    }
    return(n);
  }
  else 
    return(0);
}

/* Weighted sum of squares */

double wssq(const double *y, int n, const double *weights) {
  double res = 0.0;
  int i;
  if (weights) {
    for (i=0; i<n; i++) {
      double yi = y[i];
      res += weights[i]*yi*yi;
    }
  }
  else {
    for (i=0; i<n; i++) {
      double yi = y[i];
      res += yi*yi;
    }
  }
  return(res);
}
 
}
