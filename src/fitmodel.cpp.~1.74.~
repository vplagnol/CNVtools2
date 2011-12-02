#include "cmath"
#include "fitmodel.h"
#include "glm_test.h"
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "math.h"
#include "MyException.h"
#include "gmath.h"

extern "C" double zeroin(double ax, double bx, double (*f)(double), double tol);
static double tvalue;

vector<double> CNV_signal::GetPosterior() const 
{
  vector<double> ret;
  ret.assign(length*7, 0.);

  for (int i = 0; i != length; i++) {
    ret[ i + 0*length ] = posterior[i]; 
    ret[ i + 1*length]  = mean[i]; 
    ret[ i + 2*length]  = variance[i]; 
    ret[ i + 3*length]  = logp[i]; 
    ret[ i + 4*length]  = alpha[i]; 
    ret[ i + 5*length]  = proba_disease[i]; 
    ret[ i + 6*length]  = nu[i];
  }

  return ret;
}  

void CNV_signal::MaximizeAlpha()
{
  // Here is where the different models manifest

  if( model == HETERO ){
    if( hypothesis == H0 ) this->MaximizeAlpha(2);
    else this->MaximizeAlpha(1); 
  }
  else if( model == DISEASE ){
    this->MaximizeAlpha(2); 
    if ((nstrat_assoc > 1) || ( hypothesis == H1 )) this->MaximizeDisease();
  }
  else if( model == QT ){
    this->MaximizeAlpha(2);
    if( hypothesis == H1 ) this->MaximizeQuantitativeTrait();
  }
  else{
    Rprintf("WARNING : MaximizeAlpha not performed\n");
  }
}

void CNV_signal::MaximizeAlpha(const int& t) 
{
  // t == 1 : Separate vector of alpha for each cohort
  // t == 2 : One vector of alpha

  for (int i = 0; i != ncohorts; i++) {  //for each cohort
    for (int j = 0; j != ncomp; j++) {  //for each component
      alphas[i][j] = 0.;
    }
  }
  
  if(t==1){
    for (int i = 0; i != length; i++) {
      alphas[ cohort[i] - 1 ][ cn[i] ] += posterior[i];
    }

    vector<double> n(ncohorts,0.);
    for(int i=0; i != ncohorts; ++i){
      for (int j = 0; j != ncomp; j++){
	n[i] += alphas[i][j];
      }
    }
    for(int i=0; i != ncohorts; ++i){
      for (int j = 0; j != ncomp; j++) {
	alphas[i][j] /= n[i];
	if (alphas[i][j] < min_n/nind) alphas[i][j] = 0.;  //I want at least min_n individuals per group  
      }
    } 
    for (int i = 0; i != length; i++) alpha[i] = alphas[ cohort[i] - 1 ][ cn[i] ];
  }

  if(t==2){
    for (int i = 0; i != length; i++) {
      alphas[0][ cn[i] ] += posterior[i];
    }

    double sum = 0.; 
    for (int i = 0; i != ncomp; i++)  sum += alphas[0][i];
    for (int i = 0; i != ncomp; i++) {
      alphas[0][i] /= sum;
      if (alphas[0][i] < min_n/nind) alphas[0][i] = 0.;  //I want at least min_n individuals per group  
    }
    for (int i = 0; i != length; i++) alpha[i] = alphas[0][ cn[i] ];
  }

  
}

//glm for he probability to be a case
void CNV_signal::MaximizeDisease() 
{   
  double conv = 0.001;
  int maxit = 20;
  int init = 0;

  //// output
  double scale;
  int df_resid, rank;
  
      
  for (int i = 0; i != length; i++) {weights[i] = posterior[i];}
  
  int failure = 1;
  if (nstrat_assoc == 1) {
    
    if (hypothesis == H1) {
      failure = glm_fit(BINOMIAL, LOGIT, length, designColDisease, 1, disease_status, weights, offset, X_disease, assoc_strata, maxit, conv, init,  //input
			&rank, Xb, temp, residuals, newWeights, &scale,&df_resid);   //output
    }
    
    if (hypothesis == H0) {
      cout<<"Should not go here\n", exit(1);
      failure = glm_fit(BINOMIAL, LOGIT, length, 0, 1, disease_status, weights, NULL, X_disease, assoc_strata, maxit, conv, init,  //input
			&rank, Xb, temp, residuals, newWeights, &scale,&df_resid);   //output
    }
    
  } else {
    
    if (hypothesis == H1) {
      failure = glm_fit(BINOMIAL, LOGIT, length, designColDisease, nstrat_assoc, disease_status, weights, offset, X_disease, assoc_strata, maxit, conv, init,  //input
			&rank, Xb, temp, residuals, newWeights, &scale,&df_resid);   //output    
    }
    
    if (hypothesis == H0) {
      failure = glm_fit(BINOMIAL, LOGIT, length, 0, nstrat_assoc, disease_status, weights, NULL, X_disease, assoc_strata, maxit, conv, init,  //input. no offset under the NULL
			&rank, Xb, temp, residuals, newWeights, &scale,&df_resid);   //output    
      //for(int i = 0; i != length; i++) {if (temp[i] != temp[i]) cout<<"kkkkkkkkkkkkkkk "<<i<<"  "<<temp[i]<<"\t"<<assoc_strata[i]<<"\t"<<weights[i]<<endl;}
    }
   
  }

  if (failure == 1) {
    Rprintf("MaximizeDisease : Failure to converge\n");
  } 
  else {
    for(int i = 0; i != length; i++) proba_disease[i] = temp[i];
    this->FillGaps();
  }

}

void CNV_signal::MaximizeQuantitativeTrait()
{
  //input parameters
  int S = 1;   //only one stratum
  int maxit = 30, init = 0;  //meaningless, just OLS here
  double conv = -1.;

  //// output
  double scale;
  int df_resid, rank;

  for (int i = 0; i != length; i++) {
    weights[i] = posterior[i];
  }

  int failure = glm_fit(GAUSSIAN, IDENTITY, length, designColDisease,S, disease_status, weights, NULL, X_disease, stratum, maxit, conv, init,  //input
			&rank, Xb, temp, residuals, newWeights, &scale,&df_resid);   //output, placed in temp in that case

  if (failure == 1) {
    Rprintf("MaximizeQuantitativeTrait : Failure to converge\n");
  } 
  else {
    rms = 0;
    for(int i = 0; i != length; i++){
      proba_disease[i] = temp[i];
      rms += residuals[i]*residuals[i]*posterior[i];
    }
    rms /= nind;
    //Rprintf("Fitted rms : %4.4f \n",sqrt(rms) );

    this->FillGaps();    
  }
  
}

double CNV_signal::GetLogLikelihood () const   //meant to compute the likelihood 
{
  double logLikelihood = 0.;
    
  for (int i = 0; i != nind; i++) {   //for each individual      
    int my_max = 0;
    for (int j = 1; j < ncomp; j++) {if (logp[ j *nind + i] >= logp[ my_max *nind + i ]) {my_max = j;}}
    double tempLike = logp[ my_max *nind + i ];
      
    double temp = 0.;
    for (int j = 0; j < ncomp; j++) {temp +=  exp(logp[ j * nind + i ] - logp[ my_max*nind + i ]);}      
    logLikelihood += tempLike + log(temp);      
  }
    
  return logLikelihood;
}

//
//
//
//
/* gaussian distribution specific */
//
//
//
//

void CNV_signal::MaximizeVariancesG()
{
  //input parameters
  int S = 1;   //it means that I fit an intercept
  double conv = 0.001;
  int maxit = 40;
  int init = 0;


  //// output
  double scale;
  int df_resid, rank;

  
  for (int i = 0; i != length; i++)  { 
    weights[i] = posterior[i];
    residuals[i] = residuals[i]*residuals[i];
  }
  
  int failure;
  if (nstrat_var == 1) {
    failure = glm_fit(GAMMA, LOG, length, designColVariances, S, residuals, weights,  NULL, X_variance, stratum, maxit, conv, init,  //input
			  &rank, Xb, temp, fitted, newWeights, &scale,&df_resid);   //output
  } else {
    failure = glm_fit(GAMMA, LOG, length, 0, nstrat_var, residuals, weights,  NULL, X_variance, variance_strata, maxit, conv, init,  //input: in this version there is no design matrix, just strata
			  &rank, Xb, temp, fitted, newWeights, &scale,&df_resid);   //output
  }


  if (failure == 0) {
    for (int i = 0; i != length; i++) {variance[i] = temp[i];}
  } 
  else {
    //Rprintf("No convergence for variances\n");
  }

  FillGaps();
  for (int i = 0; i != length; i++){
    if ( (variance[i] < 0) || (variance[i] != variance[i]) ) { 
      stringstream throwmessage;
      throwmessage << "CNV_signal::MaximizeVariances : Negative variances : "
                   << i << "\t" << posterior[i] << "\t" << cohort[i] << "\t" << variance[i] << "\t" << cn[i];

      throw MyException(throwmessage.str().c_str(),1);

    }
  }
}

void CNV_signal::MaximizeMeansG()
{
  //input parameters
  int S = 1;   //only one stratum
  int maxit = 30, init = 0;  //meaningless, just OLS here
  double conv = -1.;

  //// output
  double scale;
  int df_resid, rank;

  for (int i = 0; i != length; i++)   {weights[i] = posterior[i]/variance[i];}

  int failure;
  if (nstrat_mean == 1) {
    failure = glm_fit(GAUSSIAN, IDENTITY, length, designColMeans,S, signal,weights, NULL, X_mean,stratum,maxit,conv,init,  //input
			  &rank,Xb, temp, residuals, newWeights, &scale,&df_resid);   //output, placed in "mean" in that case
  } else {
    failure = glm_fit(GAUSSIAN, IDENTITY, length, 0, nstrat_mean, signal,weights, NULL, X_mean, mean_strata,maxit,conv,init,  //input
			  &rank,Xb, temp, residuals, newWeights, &scale,&df_resid);   //output, placed in "mean" in that case    
  }

  if (failure == 0) {for (int i = 0; i != length; i++) mean[i] = temp[i];}
}

void CNV_signal::MaximizeVariancesPosteriorG(const int& t)
{
  // Note this does not yet use any linear model

  // t == 1 : fit with no diff bias
  // t == 2 : allow diff bias
 
  vector<vector<double> > n(ncohorts, vector<double>(ncomp, 0 ));

  // Fill means
  FillGaps();
  for (int i = 0; i != length; i++) means[ cohort[i] - 1][ cn[i] ] = mean[i];
  
  // Initialize
  for(int i=0; i < ncohorts; ++i){
    for(int j=0; j<ncomp; ++j){
      variances[i][j] = 0;
    }
  }

  for (int i = 0; i != length; i++){
    if(t==1){
      n[ 0 ][ cn[i] ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ 0 ][ cn[i] ] += posterior[i]*res*res;
    }
    if(t==2){
      n[ cohort[i] - 1][ cn[i] ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ cohort[i] - 1][ cn[i] ] += posterior[i]*res*res;
    }  
  }

  if(t==1){
    for(int j=0; j<ncomp; ++j){
      double temp = 0;
      for(int i=0; i < ncohorts; ++i){
	temp += shrinkage*n[ 0 ][j]*(means[i][j] - mean_p)*(means[i][j] - mean_p)/(shrinkage + n[ 0 ][j]); 
      }
      variances[ 0 ][j] = (scale + temp + variances[ 0 ][j])/( dof + n[ 0 ][j] + 3 );
    }
  }

  if(t==2){
    for(int i=0; i < ncohorts; ++i){
      for(int j=0; j<ncomp; ++j){
	double temp = scale + shrinkage*n[i][j]*(means[i][j] - mean_p)*(means[i][j] - mean_p)/(shrinkage + n[i][j]) + variances[i][j]; 
	variances[i][j] = temp/( dof + n[i][j] + 3 );
      }
    }
  }
  

  // Now fill the variance vector
  if(t==1){
    for (int i = 0; i != length; i++) variance[i] = variances[ 0 ][ cn[i] ];
  }
  if(t==2){
    for (int i = 0; i != length; i++) variance[i] = variances[ cohort[i] - 1][ cn[i] ];
  }

}

void CNV_signal::MaximizeMeansPosteriorG()
{
  // Note this does not yet use any linear model
  
  vector<vector<double> > n(ncohorts, vector<double>(ncomp, 0 ));
  
  // Initialize
  for(int i=0; i < ncohorts; ++i){
    for(int j=0; j<ncomp; ++j){
      means[i][j] = 0;
    }
  }

  for (int i = 0; i != length; i++){
    n[ cohort[i] - 1][ cn[i] ] += posterior[i];
    means[ cohort[i] - 1][ cn[i] ] += posterior[i]*signal[i];
  }

  for(int i=0; i < ncohorts; ++i){
    for(int j=0; j<ncomp; ++j){
      means[i][j] = (means[i][j] + shrinkage*mean_p)/( n[i][j] + shrinkage );
    }
  }

  // Now fill the mean vector
  for (int i = 0; i != length; i++) mean[i] = means[ cohort[i] - 1][ cn[i] ];

}

void CNV_signal::ComputePosterior ()
{

  for (int i = 0; i != nind; i++) {   //for each individual

    for (int j = 0; j != ncomp; j++) {  //for each component
      if (logp[ j*nind + i ] == -HUGE_VAL) {posterior[ j*nind + i ] = 0.;} else {
      
	double po = 1.;	
	for (int k = 0; k < ncomp; k++) {if (k != j) po += exp(logp[ k * nind + i ] - logp[ j*nind + i ]);}  //compare with others
	posterior[ j*nind + i ] = 1./po; 
	if (posterior[ j*nind + i ] < 0.0001)   {posterior[ j*nind + i ] = 0.;}
	
	if (1./po != 1./po) {
	  cerr<<po<<" makes no sense\n";
	  for (int k = 0; k < ncomp; k++) {cout<<logp[ k * nind + i ]<<endl;}
	  //for (int k = 0; k < ncomp; k++) cout<<variance[  k * nind + i ]<<"\t"<<mean[ k * nind + i ]<<"\t"<<signal[k * nind + i ]<<"\t"<<logp[i]<<"\t"<<posterior[i]<<endl;
	  exit(1);
	}
      }
    }
  }
}




void CNV_signal::ExpectationG() 
{
  double pi = 3.141592653589793238465;
  //cout<<"proba dis: "<<proba_disease[0]<<endl;
  
  for (int i = 0 ; i != length; i++) {
    double u = signal[i] - mean[i];
    double p_disease = disease_status[i] ? proba_disease[i] : 1. - proba_disease[i];
    
    if(model == DISEASE || model == HETERO){
      logp[i] = alpha[i] > 0 ? -0.5*log(2.*pi*variance[i]) - 0.5*u*u/variance[i] + log(alpha[i]) + log(p_disease) : -HUGE_VAL;
    }    
    else{
      // Gaussian contribution from QT regression
      double u_qt = disease_status[i] - proba_disease[i];
      logp[i] = alpha[i] > 0 ? -0.5*log(2.*pi*variance[i]) - 0.5*u*u/variance[i] - 0.5*log(2.*pi*rms) - 0.5*u_qt*u_qt/rms + log(alpha[i]) : -HUGE_VAL;
    }

    if (logp[i] != logp[i]) {
     
      stringstream throwmessage; 
      throwmessage << "CNV_signal::Expectation : NaN in the likelihood computation : " 
		   << i << "\t" << variance[i] << "\t" << alpha[i] << "  " << log(alpha[i]) << "  " << cn[i] << "\t" << proba_disease[i];
       
      throw MyException(throwmessage.str().c_str(),1); 
    }
  }
}

//
//
//
//
/* t distribution specific */
//
//
//
//

double CNV_signal::logpT(double x, double mu, double var, double nu)
{
  double pi = 3.141592653589793238465;
  double delta = (x-mu)*(x-mu)/var;
  double xx = (1 + nu)/2;

  double logt = gmath::lgamma(xx) - 0.5*log(var) -0.5*log(pi*nu) - gmath::lgamma(nu/2) -xx*log(1+delta/nu);
  //cout << x << "\t" << mu << "\t" << var << "\t" << nu << "\t" << logt << endl; 
  return logt;
}

void CNV_signal::ExpectationT() 
{
  double pi = 3.141592653589793238465;

  for (int i = 0 ; i != length; i++) {
    
    if(model == DISEASE || model == HETERO){
      double p_disease = disease_status[i] ? proba_disease[i] : 1. - proba_disease[i];
      logp[i] = alpha[i] > 0 ? logpT(signal[i],mean[i],variance[i],nu[i]) + log(alpha[i]) + log(p_disease) : -HUGE_VAL;
    }    
    else{
      // Gaussian contribution from QT regression
      double u_qt = disease_status[i] - proba_disease[i];
      logp[i] = alpha[i] > 0 ? logpT(signal[i],mean[i],variance[i],nu[i]) - 0.5*log(2.*pi*rms) - 0.5*u_qt*u_qt/rms + log(alpha[i]) : -HUGE_VAL;
    }

    if (logp[i] != logp[i]) {
      stringstream throwmessage; 
      throwmessage << "CNV_signal::Expectation : NaN in the likelihood computation : " 
		   << i << "\t" << variance[i] << "\t" << alpha[i] << "  " << cn[i] << "\t" << proba_disease[i];
      throw MyException(throwmessage.str().c_str(),1); 
    }

    // Need to calculate uij for t distribution
    double delta = (signal[i] - mean[i])*(signal[i] - mean[i])/variance[i];
    u[i] = (nu[i] + 1)/(nu[i] + delta);
  
  }
}

void CNV_signal::MaximizeMeansT(const int& t)
{
  // Note this does not use any linear model
  
  // t == 1 : fit with no diff bias
  // t == 2 : allow diff bias

  vector<vector<double> > n(ncohorts, vector<double>(ncomp, 0 ));
  
  // Initialize
  for(int i=0; i < ncohorts; ++i){
    for(int j=0; j<ncomp; ++j){
      means[i][j] = 0;
    }
  }


  if(t==3) {  //cn only
    
    for (int i = 0; i != length; i++){
      n[ 0 ][ cn[i] ] += posterior[i]*u[i];
      means[ 0 ][ cn[i] ] += posterior[i]*u[i]*signal[i];
    }

    for(int j=0; j<ncomp; ++j){
      if (n[0][j] > 0) means[ 0 ][j] /= n[ 0 ][j]; else means[0][j] = -99.;
    }
      
    // Now fill the mean vector
    for (int i = 0; i != length; i++) mean[i] = means[ 0 ][ cn[i] ];
    
  }

  if(t == 4) {  //with differential bias
    for (int i = 0; i != length; i++){
      n[ cohort[i] - 1][ cn[i] ] += posterior[i]*u[i];
      means[ cohort[i] - 1][ cn[i] ] += posterior[i]*u[i]*signal[i];
    }

    for(int i=0; i < ncohorts; ++i){
      for(int j=0; j<ncomp; ++j){
	if (n[i][j] > 0) means[i][j] /= n[i][j]; else means[i][j] = -99.;
      }
    }

    // Now fill the mean vector
    for (int i = 0; i != length; i++) mean[i] = means[ cohort[i] - 1][ cn[i] ];
  }


}

void CNV_signal::MaximizeVariancesT(const int& t)
{
  // Note this does not yet use any linear model

  // t == 1 : ndb const variance
  // t == 2 : db const variance
  // t == 3 : ndb free variance
  // t == 4 : db free variance
 
  vector<vector<double> > n(ncohorts, vector<double>(ncomp, 0 ));

  // Fill means
  //FillGaps();
  for (int i = 0; i != length; i++) means[ cohort[i] - 1][ cn[i] ] = mean[i];
  
  // Initialize
  for(int i=0; i < ncohorts; ++i){
    for(int j=0; j<ncomp; ++j){
      variances[i][j] = 0;
    }
  }

  if(t==1){
    for (int i = 0; i != length; i++){
      n[ 0 ][ 0 ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ 0 ][ 0 ] += posterior[i]*u[i]*res*res;
    }
    variances[ 0 ][ 0 ] = variances[ 0 ][ 0 ]/n[ 0 ][ 0 ];
    
    for (int i = 0; i != length; i++) variance[i] = variances[ 0 ][ 0 ];
  }

  if(t==2){
    for (int i = 0; i != length; i++){
      n[ cohort[i] - 1][ 0 ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ cohort[i] - 1 ][ 0 ] += posterior[i]*u[i]*res*res;
    }
    
    for(int i=0; i < ncohorts; ++i){
      if (n[i][0] >= 2) variances[ i ][ 0 ] /= n[ i ][ 0 ]; else variances[ i ][ 0 ] = 0.0001;
    }
    for (int i = 0; i != length; i++) variance[i] = variances[ cohort[i] - 1][ 0 ];
  }

  if(t==3){
    for (int i = 0; i != length; i++){
      n[ 0 ][ cn[i] ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ 0 ][ cn[i] ] += posterior[i]*u[i]*res*res;
    }
    for(int j=0; j<ncomp; ++j){
      if (n[0][j] >= 2) variances[ 0 ][j] /= n[ 0 ][j]; else variances[ 0 ][j] = 0.0001;
    }
    for (int i = 0; i != length; i++) variance[i] = variances[ 0 ][ cn[i] ];
  }

  if(t==4){
    for (int i = 0; i != length; i++){
      n[ cohort[i] - 1][ cn[i] ] += posterior[i];
      double res = signal[i] - means[ cohort[i] - 1 ][ cn[i] ];
      variances[ cohort[i] - 1][ cn[i] ] += posterior[i]*u[i]*res*res;
    }
    
    for(int i=0; i < ncohorts; ++i){
      for(int j=0; j<ncomp; ++j){
	if (n[i][j] >= 2) variances[i][j] /= n[i][j]; else variances[i][j] = 0.0001;  //this may be an issue: not enough points in a specific cohort to estimate the variance of this component
      }
    }
    for (int i = 0; i != length; i++) variance[i] = variances[ cohort[i] - 1][ cn[i] ];
  }
  return;
}

double function_to_find_root( double x ) {return -gmath::psi(0.5*x) + log(0.5*x) + tvalue;}


void CNV_signal::MaximizeNuT(const int& t)
{
  // Note this does not yet use any linear model

  // t == 1 : ndb const nu
  // t == 2 : db const nu
  // t == 3 : ndb free nu
  // t == 4 : db free nu

  vector<vector<double> > n(ncohorts, vector<double>(ncomp, 0 )); //sum(z[,i])  
  vector<vector<double> > val(ncohorts, vector<double>(ncomp, 0 )); //sum( z[,i]*(log(u[,i]) - u[,i]) )
  double numax = 25;
  double rmin = 0.1;
  double rmax = 40;

  //FillGaps();

  if(t==1){
    for (int i = 0; i != length; i++){
      n[ 0 ][ 0 ] += posterior[i];
      val[ 0 ][ 0 ] += posterior[i]*( log(u[i]) - u[i] );
    }
 
    for(int j=0; j<ncomp; ++j){
      val[ 0 ][ 0 ] = val[ 0 ][ 0 ]/n[ 0 ][ 0 ];
    }
    
    if( nus[ 0 ][ 0 ] < numax ){
      // now set tvalue = 1 + sum( z[,i]*(log(u[,i]) - u[,i]) )/sum(z[,i]) + psi(xx) - log(xx)
      double xx = (nus[ 0 ][ 0 ] + 1)/2; // uses previous value of nu
      tvalue = 1 + val[ 0 ][ 0 ] + gmath::psi(xx) - log( xx );

      int sign_start = function_to_find_root(rmin) > 0 ? +1 : -1;
      int sign_end = function_to_find_root(rmax) > 0 ? +1 : -1;

      if( sign_start != sign_end ){
	nus[ 0 ][ 0 ] =  zeroin(rmin, rmax, function_to_find_root, 0.05);
      }
      else{
	//cerr << "Warning maximizing nu: function does not cross zero" << endl;
      }

    }
    for (int i = 0; i != length; i++) nu[i] = nus[ 0 ][ 0 ];
  }

  if(t==2){
    for (int i = 0; i != length; i++){
      n[ cohort[i] - 1][ 0 ] += posterior[i];
      val[ cohort[i] - 1 ][ 0 ] += posterior[i]*( log(u[i]) - u[i] ); 
    }
    
    for(int i=0; i < ncohorts; ++i){
      val[i][ 0 ] = val[ i ][ 0 ]/n[ i ][ 0 ];
    }
  
    for(int i=0; i < ncohorts; ++i){
      if( nus[ i ][ 0 ] >= numax ) continue;

      double xx = (nus[ i ][ 0 ] + 1)/2; 
      tvalue = 1 + val[ i ][ 0 ] + gmath::psi(xx) - log( xx );

      int sign_start = function_to_find_root(rmin) > 0 ? +1 : -1;
      int sign_end = function_to_find_root(rmax) > 0 ? +1 : -1;

      if( sign_start != sign_end ){
	nus[ i ][ 0 ] =  zeroin(rmin, rmax, function_to_find_root, 0.05);
      }
      else{
	//cerr << "Warning maximizing nu: function does not cross zero" << endl;
      }

    }
    for (int i = 0; i != length; i++) nu[i] = nus[ cohort[i] - 1][ 0 ];
  }

  if(t==3){
    for (int i = 0; i != length; i++){
      n[ 0 ][ cn[i] ] += posterior[i];
      val[ 0 ][ cn[i] ] += posterior[i]*( log(u[i]) - u[i] );
    }
 
    for(int j=0; j<ncomp; ++j){
      val[ 0 ][j] = val[ 0 ][j]/n[ 0 ][j];
    }
    
    for(int j=0; j<ncomp; ++j){
      if( nus[ 0 ][ j ] >= numax ) continue;

      double xx = (nus[ 0 ][ j ] + 1)/2;
      tvalue = 1 + val[ 0 ][ j ] + gmath::psi(xx) - log( xx );

      int sign_start = function_to_find_root(rmin) > 0 ? +1 : -1;
      int sign_end = function_to_find_root(rmax) > 0 ? +1 : -1;

      if( sign_start != sign_end ){
	nus[ 0 ][ j ] =  zeroin(rmin, rmax, function_to_find_root, 0.05);
      }
      else{
	//cerr << "Warning maximizing nu: function does not cross zero" << endl;
      }
    
    }
    for (int i = 0; i != length; i++) nu[i] = nus[ 0 ][ cn[i] ];
  }

  if(t==4){
    for (int i = 0; i != length; i++){
      n[ cohort[i] - 1][ cn[i] ] += posterior[i];
      val[ cohort[i] - 1 ][ cn[i] ] += posterior[i]*( log(u[i]) - u[i] ); 
    }
    
    for(int i=0; i < ncohorts; ++i){
      for(int j=0; j<ncomp; ++j){
	if (n[i][j] > 0) val[i][j] /= n[i][j];
      }
    }
  
    for(int i=0; i < ncohorts; ++i){
      for(int j=0; j<ncomp; ++j){
	if (n[ i ][ j ] > 0) {
	  if( nus[ i ][ j ] >= numax ) continue;


	  double xx = (nus[ i ][ j ] + 1)/2;
	  tvalue = 1 + val[ i ][ j ] + gmath::psi(xx) - log( xx );
	  
	  int sign_start = function_to_find_root(rmin) > 0 ? +1 : -1;
	  int sign_end = function_to_find_root(rmax) > 0 ? +1 : -1;
	  
	  if( sign_start != sign_end ){
	    nus[ i ][ j ] =  zeroin(rmin, rmax, function_to_find_root, 0.05);
	  }
	  else{
	    //cerr << "Warning maximizing nu: function does not cross zero" << endl;
	  }
	} else { nus[ i ][ j ] = 10;}  //some default values when proba of class is 0
      }
    }
    for (int i = 0; i != length; i++) nu[i] = nus[ cohort[i] - 1][ cn[i] ];  
  }

  return;
}


void CNV_signal::PrintOneLine(const int i) const 
{
  cout<<i<<"\t"<<signal[i]<<"\t"<<mean[i]<<"\t"<<residuals[i]<<"\t"<<posterior[i]<<"\t"<<variance[i]<<"\t"<<cn[i]<<"\t"<<logp[i]<<endl;
}

void CNV_signal::Print() const 
{        
  for (int i = 0; i != ncomp; i++) {
    cout<<"Component "<<i<<endl;
    cout<<"Mean: "<<mean[0 + i*nind]<<"\tStd. dev:"<<variance[0 + i*nind]<<"\talpha:"<<alpha[0 + i*nind]  <<endl;
  }

  //for (int i = 0; i != ncomp; i++) {cout<<mean[0 + i*nind]<<"\t"<<variance[0 + i*nind]<<"\t"<<logp[0+i*nind]<<"\t"<<signal[0+nind*i]<<"\t"<<posterior[ 0 + nind*i]<<endl;}
  cout<<"\n\n\n";
    
}

void CNV_signal::PrintParams() const
{

  for(int j=0; j<ncomp; ++j){
    for(int i=0; i < ncohorts; ++i){
      cout << "\t" << nus[i][j]; 
    }
    //cout << "\n";
  }
  cout << "\n";
}

void CNV_signal::Check_order()
{
  vector<double> alpha_m(ncomp, 0);

  for(int i=0; i != ncohorts; ++i){

    for(int j=0; j != ncomp; ++j){
      nus[i][j] = -99.0;
      variances[i][j] = -99.;
      means[i][j] = -99.;
      postLogit[j] = 0.;
      postTable[i][j] = -1.;
    }
  }

  
  for(int s = 0; s != nstrat_assoc; ++s){
    for(int j=0; j != ncomp; ++j){
      postLogit2[ s ][ j ] = -0.1;
      proba_d[ s ][ j ] = 0.;
    }
  }



  for (int i = 0; i != length; i++) {
    
    if (posterior[i] > postTable [cohort[i] - 1 ][ cn[i] ] ) {
      variances[ cohort[i] - 1 ][ cn[i] ] = variance[i];
      means[ cohort[i] - 1 ][ cn[i] ] = mean[i];
      postTable [cohort[i] - 1 ][ cn[i] ] = posterior[i];
      nus[ cohort[i] - 1 ][ cn[i] ] = nu[i];
    }
   
    int s =  assoc_strata[i] - 1;
    if (posterior[i] > postLogit2[ s ][ cn[i] ]) {
      proba_d[ s ][ cn[i] ] = proba_disease[i];
      postLogit2[ s ][ cn[i] ] = posterior[i];
    }

    if (posterior[i] > postLogit[ cn[i] ]) {
      postLogit[ cn[i] ] = posterior[i];
      alpha_m[ cn[i] ] = alpha[i];
    }
  }

  vector<int> tmp (ncomp, -99);
  vector<vector<int> > ranks;
  for(int i=0; i != ncohorts; ++i){
    myRank * hope = new myRank(means[ i ]);
    hope->get_orders(tmp);
    ranks.push_back(tmp);
    //for (int j = 0; j != ncomp; j++) cout<<ranks[i][j]<<"\t";cout<<endl;
    delete hope;
  }

  for(int i=0; i != ncohorts; ++i){
    for (int j = 0; j != length; j++) {
      variance[j] = variances[ cohort[j] - 1 ][ ranks[cohort[j] - 1][ cn[j] ] ];
      nu[j] = nus[ cohort[j] - 1 ][ ranks[cohort[j] - 1][ cn[j] ] ];
      mean[j] = means[ cohort[j] - 1 ][    ranks[cohort[j]-1][cn[j]] ];
      alpha[j] = alpha_m [  ranks[cohort[j]-1][cn[j]]   ];

      proba_disease[j] = proba_d[ assoc_strata[j] - 1 ][    ranks[cohort[j]-1][cn[j]]  ];
    }
  }
} 

void CNV_signal::FillGaps() 
{
  
  for(int i = 0; i != ncohorts; ++i){
    for(int j = 0; j != ncomp; ++j){
      variances[i][j] = 0.0001;
      nus[i][j] = -99.;
      means[i][j] = -99.;
      postLogit[j] = 0.;
      postTable[i][j] = 0;
    }
  }

  for(int i = 0; i != nstrat_assoc; ++i){
    for(int j = 0; j != ncomp; ++j){
      postLogit2[ i ][ j ]  = -0.1;
    }
  }


  for (int i = 0; i != length; i++) {
    if (posterior[i] > postTable [cohort[i] - 1 ][ cn[i] ] ) {
      nus[ cohort[i] - 1 ][ cn[i] ] = nu[i];
      if (variance[i] > 0) variances[ cohort[i] - 1 ][ cn[i] ] = variance[i];
      means[ cohort[i] - 1 ][ cn[i] ] = mean[i];
      postTable [cohort[i] - 1 ][ cn[i] ] = posterior[i];
    }
   
    if (posterior[i] > postLogit[ cn[i] ]) {
      postLogit[ cn[i] ] = posterior[i];
    }
    
    int s = assoc_strata[ i ] - 1;
    if (posterior[i] > postLogit2[ s ][ cn[i] ]) {
      if (proba_disease[i]  != proba_disease[i]) cout<<"Bug when computing disease probability "<<proba_disease[i]<<"\t"<<posterior[i]<<"\t"<<cn[i]<<"\t"<<s<<"\t"<<hypothesis<<endl;
      proba_d[ s ] [ cn[i] ] = proba_disease[i];
      postLogit2[ s ][ cn[i] ] = posterior[i];      
    }

  }

  //PrintParams();

  for (int i = 0; i != length; i++) {   //now fill the situations where the weight is 0 with the situations where the weight is not 0!!
    nu[i] = nus[ cohort[i] - 1 ][ cn[i] ];
    variance[i] = variances[ cohort[i] - 1 ][ cn[i] ];
    mean[i] = means[ cohort[i] - 1 ][ cn[i] ];
    proba_disease[i] = proba_d[ assoc_strata[i] - 1 ][ cn[i] ];
    //if (proba_disease[i]  != proba_disease[i]) cout<<"llllll "<<proba_disease[i]<<"\t"<< proba_d[ assoc_strata[i] - 1 ][ cn[i] ]<<endl;
  }

}
  

CNV_signal::CNV_signal(const int nind_a, 
		       const int ncomp_a, 
		       const int * cohort_a, 
		       const double * signal_a, 
		       const double * disease, 
		       const double * mean_start_a, 
		       const double * var_start_a, 
		       const double * nu_start_a, 
		       const double * alpha_start_a,
		       const double * offset_a,
		       const double * mean_design, 
		       const double * var_design, 
		       const double * disease_design, 
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
		       ) : nind (nind_a), 
			   ncomp(ncomp_a), 
			   cohort(cohort_a),
			   signal(signal_a), 
			   disease_status(disease), 
			   offset(offset_a),
			   X_mean(mean_design), 
			   X_variance(var_design), 
			   X_disease(disease_design), 
			   designColMeans(designColMeans_a), 
			   designColVariances(designColVariances_a), 
			   designColDisease(designColDisease_a),
			   model(m),
			   hypothesis(h),
			   min_n(min_n_a),
			   variance_strata ( variance_strata_a ),
			   nstrat_var (nstrat_var_a),
			   mean_strata ( mean_strata_a ),
			   nstrat_mean (nstrat_mean_a),
			   assoc_strata ( assoc_strata_a ),
			   nstrat_assoc (nstrat_assoc_a)
			   
{    
  length = nind*ncomp;    
  ncohorts = *max_element(cohort, cohort + length);
    
  ///////// now make space for the arrays that will contain the parameters
  vector<double> def(ncomp, 0.0);
  postLogit.assign(ncomp, 0);

  for (int i = 0; i != nstrat_assoc; i++) {
    proba_d.push_back(def);
    postLogit2.push_back(def);
  }

  for (int i = 0; i != ncohorts; i++) {
    variances.push_back(def);
    means.push_back(def);
    postTable.push_back(def);
    alphas.push_back(def);
    nus.push_back(def);
  }
  
  temp = new double [length];
  mean = new double [length];
  variance = new double [length];
  nu = new double [length];
  alpha = new double [length];
  logp = new double [length];
  posterior = new double [length];
  u = new double[length];
  fitted = new double [length];
  residuals = new double [length];
  weights = new double [length];
  newWeights = new double [length];
  proba_disease = new double [length];
  Xb = new double [length*max(max(designColMeans, designColVariances),designColDisease) ];        
  stratum = new int [length];

  int count = 0;
  cn = new int [length];
  for (int i = 0; i != ncomp; i++) {
    for (int j = 0; j != nind; j++) {
      cn[count] = i;
      count++;
    }
  }

  double c_disease = 0.;
  for (int i = 0; i != length; i++) {
    posterior[i] = 1./ncomp;   //initialize with non-zero values, otherwise there will be a bug
    u[i] = 1;
    stratum[i] = 1;
    mean[i] = mean_start_a[i]; 
    variance[i] = var_start_a[i];
    nu[i] = nu_start_a[i];
    c_disease += (disease_status[i] == 1);
    alpha[i] = alpha_start_a[i]; 
  }

  if(model == DISEASE || model == HETERO){
    double proba_disease_null = c_disease/length;   //probability to be a case under the null
    for (int i = 0; i != length; i++) {proba_disease[i] = proba_disease_null;}
  }
  else if(model == QT){

    //calculate mean and sigma of QT
    double mean_qt = 0;
    double sigma_qt = 0;
    
    for(int i=0; i<length; ++i){
      if(cn[i] == 0){
	mean_qt += disease_status[i];
	sigma_qt += disease_status[i]*disease_status[i];
      }
    }
  
    mean_qt /= nind;
    sigma_qt = sigma_qt/nind - mean_qt*mean_qt; 
    
    for (int i = 0; i != length; i++) {
      proba_disease[i] = mean_qt;
    }
    rms = sigma_qt;

    //Rprintf("QT params under H0 : %4.4f \t %4.4f \n",mean_qt,sigma_qt);
  }

  // Fill the prior variables
  // First calculate mean and variance
  mean_p = 0;
  double dvar = 0;
  for(int i=0; i<length; ++i){
    if(cn[i] == 0){
      mean_p += signal[i];
      dvar += signal[i]*signal[i];
    }
  }

  mean_p /= nind;
  dvar = dvar/nind - mean_p*mean_p; 
  scale = dvar/(ncomp*ncomp);
  shrinkage = 0.01;
  dof = 3;

  //Rprintf("\n**** Instantiated new CNV model ***\n");
  //Rprintf("N components   : %d \n", ncomp);
  //Rprintf("N individuals  : %d \n", nind);
  //Rprintf("N cohorts      : %d \n", ncohorts);

  //Rprintf("\n**** Fit parameters ***************\n");
  //Rprintf("col mean : %d , col var : %d \n\n\n", designColMeans, designColVariances );
  
  //Rprintf("Mean design matrix\n");
  //for(int i=0; i<length; ++i){
  //Rprintf("%d",i);
  //for(int j=0; j<designColMeans; ++j){
  //  Rprintf("\t %4.2f",X_mean[i + j*length]);
  //}
  //Rprintf("\n");
  //}

  //Rprintf("\n**** Prior parameters ************\n");
  //Rprintf("data mean : %4.4f \ndata var : %4.4f\n",mean_p,dvar);
  //Rprintf(" mean_p : %4.4f \n scale : %4.4f \n shrinkage : %4.4f \n dof : %4.2f \n\n\n", mean_p, scale, shrinkage, dof);

}


CNV_signal::~CNV_signal() 
{
  delete [] temp;
  delete [] mean;
  delete [] variance;
  delete [] nu;
  delete [] alpha;
  delete [] logp;
  delete [] fitted;
  delete [] residuals;
  delete [] weights;
  delete [] newWeights;
  delete [] Xb;
  delete [] stratum;
  delete [] cn;
  delete [] proba_disease;
  delete [] posterior;
  delete [] u;
}


