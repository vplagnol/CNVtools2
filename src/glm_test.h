#ifndef _GLM_H_
#define _GLM_H_

/* Family */

extern "C" {
#define BINOMIAL  1
#define POISSON   2
#define GAUSSIAN  3
#define GAMMA     4

/* Link */

#define LOGIT     1
#define LOG       2
#define IDENTITY  3
#define INVERSE   4

/* GLM definition functions */
  int wcenter(const double *y, int n, const double *weight, const int *stratum, int nstrata, int resid, double *ynew);
  int wresid(const double *y, int n, const double *weight, const double *x, double *ynew);
  
  double wssq(const double *y, int n, const double *weight);
  double varfun(int, double);
  int muvalid(int, double);
  double linkfun(int, double);
  double invlink(int, double);
  double dlink(int, double);
  
  

/* Fit a base model */

  int glm_fit(int family, int link, int N, int M, int S,
	      const double *y, const double *prior, const double * offset, const double *X, 
	      const int *stratum, int maxit, double conv, int init, 
	      int *rank, double *Xb, 
	      double *fitted, double *resid, double *weights, 
	      double *scale, int *df_resid);
  
  /* Score test for additional terms */
  
  void glm_score_test(int N, int M, int S, const int *stratum, 
		      int P, const double *Z, int C, const int *cluster,
		      const double *resid, const double *weights, 
		      const double *Xb, double scale,
		      double max_r2, double *chi2, int *df);

  
}
#endif
