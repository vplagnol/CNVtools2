#include "cmath"
#include "fitmodel.h"
#include <cstdlib>
#include <R.h>
extern "C" {

#include "glm_test.h"



/* Fit GLM, possibly with a stratification in the RHS

Input:

family       GLM family (see below)
link         Link function (see below)
N            # units
M            # X variables
S            # strata (0 means no intercept)
y            y-variable (N-vector)
prior        prior weights (if present)
X            If M>0, N*M matrix of X variables
stratum      If S>1, stratum assignments coded 1...S (N-vector)
maxit        Maximum number of iterations of IRLS algorithm
conv         Proportional change in weighted sum of squares residuals to
             declare convergence
init         If true (non-zero), the iteration starts from initial estimates 
             of fitted values (see below). This option has no effect if
	     no iteration is required

Output:

rank         rank of X after regression on strata
Xb           orthogonal basis for X space (N*rank matrix)
fitted       fitted values 
resid        working residuals (on linear predictor scale) (N-vector)
weights      weights (N-vector)
scale        scale factor (scalar)
df_resid     residual degrees of freedom

Return

0            convergence
1            no convergence after maxit iterations

*/

int glm_fit(int family, int link, int N, int M, int S,
	    const double *y, const double *prior, const double * offset, const double *X, 
	    const int *stratum, int maxit, double conv, int init, 
	    int *rank, double *Xb, 
	    double *fitted, double *resid, double *weights, 
	    double *scale, int *df_resid) {
  const double eta = 1.e-8;       /* Singularity threshold */
  int i = 0, j=0;
  int Nu, dfr, irls;
  int empty = 0;

  /* Is iteration necessary? */


  
  irls =  ( ((offset) || (M>0)) && !((family==GAUSSIAN) && (link==IDENTITY)));
  //if (family == BINOMIAL) cout<<"M  "<<M<<endl;
  //for (int i = 0; i != 1000; i++) {cout<<i<<"\t"<<prior[i]<<endl;}

  if (!init || !irls) {
    /* Fit intercept and/or strata part of model */
    empty = wcenter(y, N, prior, stratum, S, 0, fitted);    
  }

	

  Nu = 0;
  int invalid = 0;
  for (i=0; i<N; i++) {
    double mu = fitted[i];
    double ri, wi;
    double pi = prior? prior[i] : 1.0;
     if (!muvalid(family, mu)) {
       invalid = 1;
       pi = 0.0;
    }
    if (!(pi)) {wi = ri = 0.0;}
    else {
      Nu ++;
      double Vmu = varfun(family, mu);
      if (link == family) {
	ri = (y[i] - mu)/Vmu;
	wi = pi*Vmu;
      }
      else {
	double D = dlink(link, mu);
	ri = D*(y[i] - mu);
	wi = pi/(D*D*Vmu);
      }
    }
    weights[i] = wi;
    resid[i] = ri;
    if (weights[i] < 0.0001) weights[i] = 0.;
  }


  /* If M>0, include covariates */
  int x_rank = 0, convg = 0, iter = 0;
  if ((M == 0) && !offset) convg = 1;
  if (M> 0 || offset) {   //maybe also where there is an offset?
    convg = 0;
    double wss_last = 0.0;
    if (irls) {
      
      /* IRLS algorithm */
      double *yw = (double *) Calloc(N, double);  //working y
      while(iter<maxit && !convg) {
	for (i=0; i<N; i++) {
	  yw[i] = resid[i] + linkfun(link, fitted[i]);   //current estimate of eta + (y-mu)/gradient
	}

	if (offset) {for (i=0; i<N; i++) {yw[i] -= offset[i];}}
	empty = wcenter(yw, N, weights, stratum, S, 1, resid);    //removes the mean from yw

	//////////// now it tries to fit the regression line (no intercept) to the residuals
	const double *xi = X;
	double *xbi = Xb;
	x_rank = 0;

	for (i=0; i<M; i++, xi+=N) {
	  double ssx = wssq(xi, N, weights);
	  wcenter(xi, N, weights, stratum, S, 1, xbi);
	  double *xbj = Xb;
	  for (j=0; j<x_rank; j++, xbj+=N) wresid(xbi, N, weights, xbj, xbi);
	  double ssr = wssq(xbi, N, weights);
	  if (ssr/ssx > eta) {
	    wresid(resid, N, weights, xbi, resid);   //takes the residuals after fitting the regression line (no intercept) to the mean value per stratum
	    x_rank++;
	    xbi+=N;
	  }
	}


	double wss=0.0;
	Nu = 0;
	for (i=0; i<N; i++) {
	  double D, Vmu, ri, wi;
	  double mu = invlink(link, yw[i] - resid[i]);   //ie. (yw - (yw - mean(yw))) = mean(yw)
	  if (offset) {mu = invlink(link, yw[i] + offset[i] - resid[i]);}

	  fitted[i] = mu;

	  double pi = prior? prior[i] : 1.0;
	  if (!(pi && (weights[i]>0.0))) {wi = ri = 0.0;} else {
	    
	    if (!(muvalid(family, mu))) {
	      if ((family == 4) && (mu > 5.0)) {mu = fitted[i] = 5.0;}
	      if ((family == 4) && (mu < 0.001)) {mu = fitted[i] = 0.001;}
	    }
	    
	    Vmu = varfun(family, mu);
	    Nu ++;
	    if (link == family) {
	      ri = (y[i] - mu)/Vmu;
	      wi = pi*Vmu;
	    }
	    else {
	      D = dlink(link, mu);
	      ri = D*(y[i] - mu);
	      wi = pi/(D*D*Vmu);
	    }
	    wss += wi*ri*ri;
	    
	      
	    weights[i] = wi;
	    resid[i] = ri;
	    if (weights[i] < 0.0001) weights[i] = 0.;
	  }
	}
	    
	convg = (family==2) || (Nu<=0) || (iter && (fabs(wss-wss_last)/wss_last < conv));
	wss_last = wss;
	iter ++;
      }
      Free(yw);
    } else {  

      /* Simple linear Gaussian case */

      const double *xi = X;
      double *xbi = Xb;
      x_rank = 0;
      for (i=0; i<M; i++, xi+=N) {
	double ssx = wssq(xi, N, weights);
	wcenter(xi, N, weights, stratum, S, 1, xbi);
	double *xbj = Xb;
	for (j=0; j<x_rank; j++, xbj+=N)  wresid(xbi, N, weights, xbj, xbi);
	double ssr = wssq(xbi, N, weights);
	if (ssr/ssx > eta) {
	  wresid(resid, N, weights, xbi, resid);
	  x_rank++;
	  xbi+=N;
	}
      }
      for (i=0; i<N; i++) {fitted[i] = y[i] - resid[i];}   //need to obtain the fitted values
      wss_last = wssq(resid, N, weights);
    }
    dfr = Nu  - S + empty - x_rank;
    if (family>2) 
      *scale = wss_last/(dfr);
    else
      *scale = 1.0;
  }
  else {
    if ((S>1) && invalid) { /* Need to recalculate empty stratum count  */
      empty = wcenter(fitted, N, weights, stratum, S, 0, fitted); 
    }
    dfr = Nu - S + empty;
    if (family>2) 
      *scale = wssq(resid, N, weights)/(dfr);
    else
      *scale = 1.0;    
    x_rank = 0;
  }
  *df_resid = dfr>0? dfr : 0;
  *rank = x_rank;

  return(irls && !convg);
}


/* 

Variance function

family:
  1    Binomial
  2    Poisson
  3    Gaussian
  4    Gamma

*/
      
double varfun(int family, double mu){
  switch (family) {
  case 1: return((mu*(1.0-mu)));  /* Binomial */
  case 2: return(mu);             /* Poisson */
  case 3: return(1.0);            /* Gaussian */
  case 4: return(mu*mu);          /* Gamma */
  default: return(0.0);
  }
}

/* Valid values for fitted value, mu. 

If, during iteration, an invalid value is returned, the case is omitted 

*/

int muvalid(int family, double mu) {
  const double minb = 0.0001, maxb = 0.9999, minp = 0.0001, gammaMax = 5., gammaMin = 0.001;
  switch (family) {
  case 1: return(mu>minb && mu<maxb);    /* Binomial */
  case 2: return(mu>minp);               /* Poisson */
  case 4: return ( (mu> gammaMin) && (mu < gammaMax));               /* Gamma */
  default: return(1);                     
  }
}

/* Link function

Link
  1    Logit
  2    Log
  3    Identity
  4    Inverse

Note that a canonical link shares the code of the corresponding family
so that the test for canonical link is (link==family)

*/


double linkfun(int link, double mu) {
  switch (link) {
  case 1: {     /* Logit */
    if (mu == 1) return HUGE_VAL;
    if (mu == 0) return -HUGE_VAL;
    return(log(mu/(1.0-mu)));
  }
  case 2: return(log(mu));             /* Log */
  case 3: return(mu);                  /* Identity */
  case 4: return(-1.0/mu);              /* Inverse */
  default: return 0.0;
  }
}

double invlink(int link, double eta) {
  switch (link) {
  case 1: {  /* Logit */
    if (eta == HUGE_VAL) return (1); 
    if (eta == -HUGE_VAL) return (0); 
    return(exp(eta)/(1.0+exp(eta))); 
  }
  case 2: return(exp(eta));                /* Log */
  case 3: return(eta);                     /* Identity */
  case 4: return(-1.0/eta);                 /* Inverse */
  default: return(0.0);
  }
}

double dlink(int link, double mu) {
  switch (link) {
  case 1: return(1.0/(mu*(1.0-mu)));     /* Logit */
  case 2: return(1.0/mu);                /* Log */
  case 3: return(1.0);                   /* Identity */
  case 4: return(1.0/(mu*mu));           /* Inverse */
  default: return 0.0;
  }
}


/* 

GLM score test 

Input:

P         Number of new explanatory variables to be added 
Z         N*P matrix containing these (destroyed)
C         If robust variance estimate to be used, number of clusters
          (if C==1, each unit forms a cluster)
cluster   If C>1, cluster assignments code 1...C (N-vector)
max_R2    For P>1, the maximum value of R^2 between each column and revious 
          columns (after regression on X and strata)

For all other input arguments, see glm_fit, but note that M now coincides 
with rank -- the number of columns in Xb

Output:

chi2  Score test 
df    Degrees of freedom for asymptotic chi-squared distribution


*/

void glm_score_test(int N, int M, int S, const int *stratum, 
		    int P, const double *Z, int C, const int *cluster,
		    const double *resid, const double *weights, 
		    const double *Xb, double scale,
		    double max_R2, double *chi2, int *df) {
  const double eta1 = 1.e-8;   /* First stage singularity test */
  double eta2 = 1.0 - max_R2;   /* Second stage singularity test */
  const double *Zi = Z;
  double *Zr, *Zri;
  double *U = NULL, *Ui = NULL;
  double test = 0.0;
  int i = 0, j = 0;

  /* Work array */

  Zri = Zr = (double *) Calloc(N*P, double);
  int nc = 0;
  if (C) {
    nc = (C==1)? N: C;
    Ui = U = (double *) Calloc(nc*P, double);
  }
    

  /* Main algorithm */
 
  int rank = 0;
  for (i=0; i<P; i++, Zi+=N) {
    /* Regress each column of Z on strata indicators and X basis */
    double ssz = wssq(Zi, N, weights);
    wcenter(Zi, N, weights, stratum, S, 1, Zri);
    const double *Xbj = Xb;
    for (j=0; j<M; j++, Xbj+=N) 
      wresid(Zri, N, weights, Xbj, Zri);
    double ssr = wssq(Zri, N, weights);
    if (ssr/ssz > eta1) {     /* First singularity test */
      if (C) {
        /* Add new column to matrix of score contributions */
	if (C==1) {
	  for (j=0; j<N; j++)
	    Ui[j] = Zri[j]*weights[j]*resid[j];
	}
	else {
	  for (j=0; j<nc; j++)
	    Ui[j] = 0.0;
	  for (j=0; j<N; j++) {
	    int ic = cluster[j] - 1;
	    Ui[ic] += Zri[j]*weights[j]*resid[j];
	  }
	}
	/* Regress on previous columns */
	ssr = wssq(Ui, nc, NULL);
	double *Uj = U;
	for (j=0; j<rank; j++, Uj+=nc) 
	  wresid(Ui, nc, NULL, Uj, Ui);
	/* Sum and sum of squares */
	double sU = 0.0, ssU = 0.0;
	for (j=0; j<nc; j++) {
	  double Uij = Ui[j];
	  sU += Uij;
	  ssU += Uij*Uij;
	}
	/* Second singularity test */
	if (ssU/ssr > eta2) {
	  test += sU*sU/ssU;
	  rank++;
	  Zri += N;
	  Ui += nc;
	}
      }
      else {
	double *Zrj = Zr;
	for (j=0; j<rank; j++, Zrj+=N)
	  wresid(Zri, N, weights, Zrj, Zri);
	/* Sums and sums of squares */
	double ws = 0.0, wss = 0.0;
	for (j=0; j<N; j++) {
	  double Zrij = Zri[j];
	  double wz = weights[j]*Zrij;
	  ws += wz*resid[j];
	  wss += Zrij*wz;
	}
	/* Second singularity test */
	if (wss/ssr > eta2) {
	  test += ws*ws/(scale*wss);
	  rank++;
	  Zri += N;
	}
      }
    }
  }
  *chi2 = test;
  *df = rank;
  Free(Zr);
  if (C)
    Free(U);
}
}
