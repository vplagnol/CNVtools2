#ifndef _FITMODEL_H_
#define _FITMODEL_H_


#include<iostream>
#include<vector>
#include<string>
#include <algorithm>

#include <Rinternals.h>

enum MODEL{ DISEASE=0, HETERO=1, QT=2 };
enum HYPOTHESIS{ H0, H1 };

using namespace std;

extern "C"{
  
class CNV_signal
{
 public:
  int nind, ncomp, length, ncohorts;
  double likelihood, rms;


  double * temp;

  
  double * posterior;
  double * u;
  double * logp;   //loglikelihood
  double * weights;
  double * newWeights;
  double * Xb;
  double * proba_disease;
  int * cn;
  int * individual;

  // Prior values
  double mean_p;
  double shrinkage;
  double dof;
  double scale;

  double * fitted;
  double * residuals;
  int * stratum;

    
  const int * cohort;
  const double * signal;
  const double * disease_status;
  double * mean;
  double * variance;
  double * nu;
  double * alpha;
  const double * offset;


  const double * X_mean;
  const double * X_variance;
  const double * X_disease;

  int designColMeans, designColVariances, designColDisease;
  MODEL model;
  HYPOTHESIS hypothesis;
  double min_n;
  const int * variance_strata;
  int nstrat_var;

  const int * mean_strata;
  int nstrat_mean;

  const int * assoc_strata;
  int nstrat_assoc;





  vector<double> max_logP, proba_not_outlier; //for each individual stores the maximum logp
  double logP_threshold;
  //vector<double> postLogit;
  //vector< vector<double> >  postLogit2, proba_d;
  
  vector< vector<double> > variances, means, alphas, nus, postTable;   //depends on batch and copy number





  

  CNV_signal (const int nind, 
	      const int ncomp, 
	      const int * cohort, 
	      const double * signal, 
	      const double * disease, 
	      const double * mean, 
	      const double * variance, 
	      const double * nu,
	      const double * alpha, 
	      const double * offset,
	      const double * design_mean, 
	      const double * design_var, 
	      const double * design_disease, 
	      const int designColMeans_a, 
	      const int designColVariances_a, 
	      const int designColDisease_a,
	      MODEL m, 
	      HYPOTHESIS h, 
	      const double min_n_a,
	      const int * variance_strata_a,
	      const int nstrat_var_a,
	      const int * mean_strata_a,
	      const int nstrat_mean_a,
	      const int * assoc_strata_a,
	      const int nstrat_assoc_a
	      );

  ~CNV_signal ();

  double GetLogLikelihood() const;    
  vector<double> GetPosterior() const;
  void ComputePosterior();

  // alpha/disease model
  void MaximizeAlpha(const int& t);
  void MaximizeAlpha();
  void MaximizeDisease();
  void MaximizeQuantitativeTrait();


  // Gaussian specific functions
  void ExpectationG();
  void MaximizeMeansG();
  void MaximizeVariancesG();
  void MaximizeMeansPosteriorG();
  void MaximizeVariancesPosteriorG(const int& t);

  // t distribution specific functions
  double logpT(double x, double mu, double var, double nu);
  void ExpectationT();
  void MaximizeMeansT(const int& t);
  void MaximizeVariancesT(const int& t);
  void MaximizeNuT(const int& t);

  void FillGaps();
  void Check_order();
  
  void Print() const;
  void PrintOneLine(const int i) const;
  void PrintParams() const;
};

}

class myRank {
 public:
  vector<double> index;
  
  int operator()(int i1, int i2) const { 
    return(index[i1] < index[i2]);
  }
  
  void get_orders(vector<int> & w) const {
    w.assign(index.size(), 0);
    for (int i = 0; i != (int) index.size(); i++) w[i] = i;
    std::sort(w.begin(), w.end(), *this);
  }
  myRank (vector<double> & index_a): index (index_a) {};
};

#endif

