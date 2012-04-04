#include "cmath"
#include "fitmodel.h"
#include "MyException.h"
#include <cstdlib>

extern "C" {
  SEXP C_fitmodel (const SEXP ncomp_a, const SEXP nind_a, const SEXP hyp, const SEXP data, const SEXP offset_logit,
		   const SEXP design_mean, const SEXP design_variance, const SEXP design_disease, 
		   const SEXP control_parameters, const SEXP mix_model, const SEXP pi_model);
  
  SEXP getListElement(SEXP list, const char *str);
  SEXP get_first_match ( const SEXP length, const SEXP col1, const SEXP col2, const SEXP val1, const SEXP val2, const SEXP values);
}

/* 
   List of models available
   Gaussians
   10 : fit gaussian mixture, use design matrices to specify models
   11 : fit gaussian mixture bayes regularization without diff bias in variance
   12 : fit gaussian mixture bayes regularization wit diff bias in variance
 
   T  mean model    var model           nu model
   2     1 (ndb)       1  (ndb const)      1  (ndb const) 
         2 (db)        2  (db  const)      2  (db  const)
	               3  (ndb free)       3  (ndb free)
		       4  (db  free)       4  (db  free)

   eg 2111 : fit t, no diff bias mean, const var, const nu
*/


void fit_model_gaussian(CNV_signal* model, vector<double>& posterior, string& status, const double& max_iter, const double& tol, int& model_type);
void fit_model_t(CNV_signal* model, vector<double>& posterior, string& status, const double& max_iter, const double& tol, int& model_type);

SEXP get_first_match ( const SEXP length, const SEXP col1, const SEXP col2, const SEXP val1, const SEXP val2, const SEXP values) {
  

  const double * c1 = REAL(col1);
  const double * c2 = REAL(col2);
  const double v1 = *REAL(val1);
  const double v2 = *REAL(val2);
  const double * val = REAL(values);
  const int le = *INTEGER(length);
  
  for (int i = 0; i != le; i++) {
    if ( (c1[i] == v1) && (c2[i] == v2) )  return (ScalarReal(val[i]));    
  }
  
  return (ScalarReal(-99.));
}


SEXP C_fitmodel (const SEXP ncomp_a, 
		 const SEXP nind_a, 
		 const SEXP hyp, 
		 const SEXP data, 
		 const SEXP offset_logit,
		 const SEXP design_mean, 
		 const SEXP design_variance, 
		 const SEXP design_disease,
		 const SEXP control_parameters, 
		 const SEXP mix_model, 
		 const SEXP pi_model)
{  
  if (TYPEOF(hyp) != STRSXP) {cerr<<"hyp should be a character\n";exit(1);}
  if (TYPEOF(ncomp_a) != INTSXP)  {cerr<<"Argument error - ncomp"<<endl;exit(1);}
  if (TYPEOF(nind_a) != INTSXP)  {cerr<<"Argument error - nind"<<endl;exit(1);}
  if (TYPEOF(data) != VECSXP)  {cerr<<"Argument error - data"<<endl;exit(1);}
  if (TYPEOF(control_parameters) != VECSXP)  {cerr<<"Argument error - control_parameters"<<endl;exit(1);}
  if (TYPEOF(mix_model) != INTSXP)  {cerr<<"Argument error - mix.model"<<endl;exit(1);}
  if (TYPEOF(pi_model) != INTSXP)  {cerr<<"Argument error - pi.pmodel"<<endl;exit(1);}

  const double tol = *REAL(getListElement(control_parameters, "tol"));
  const double max_iter = *REAL(getListElement(control_parameters, "max.iter"));
  const double min_freq = *REAL(getListElement(control_parameters, "min.freq"));
  const double logP_threshold = *REAL(getListElement(control_parameters, "logP.outliers"));

  const double * logit_offset_p = REAL(offset_logit);
  const double * mean_design = REAL(design_mean);
  const double * var_design = REAL(design_variance);
  const double * disease_design = REAL(design_disease);
  const int * mean_dims = INTEGER(getAttrib(design_mean, R_DimSymbol));
  const int * variance_dims = INTEGER(getAttrib(design_variance, R_DimSymbol));
  const int * disease_dims = INTEGER(getAttrib(design_disease, R_DimSymbol));

  //cout<<"disease dims "<<disease_dims[0]<<" "<<disease_dims[1]<<endl;

  const string test (CHAR ( STRING_ELT (hyp, 0)));

  ///////////// Now I believe that all arguments are read only    
  const int ncomp = *INTEGER(ncomp_a);
  const int nind  = *INTEGER(nind_a);

  
  //const int * strat_assoc = INTEGER(getListElement(data, "strata.association"));
  const int * strat_assoc = NULL;
  const int * strat_var = INTEGER(getListElement(data, "strata.var"));
  const int * strat_mean = INTEGER(getListElement(data, "strata.mean"));
  const int * cohort = INTEGER(getListElement(data, "batch"));

 
  const double * alpha_start = REAL(getListElement(data, "alpha.start"));
  const double * disease = REAL(getListElement(data, "trait"));
  const double * signal = REAL(getListElement(data, "signal"));
  const double * nu_start = REAL(getListElement(data, "nu.start"));
  const double * mean_start = REAL(getListElement(data, "mean.start"));
  const double * var_start = REAL(getListElement(data, "var.start"));




  
  int nstrat_var = 0;
  vector<int> array_strat_var (500, 0);
  for (int i = 0; i != nind*ncomp; i++) {
    if (strat_var[i] > 499) {cerr<<"No more than 500 strata are allowed\n";exit(1);}
    array_strat_var[ strat_var[i] ]++;
    if ( array_strat_var[ strat_var[i] ] == 1 ) nstrat_var++;
  }

  int nstrat_mean = 0;
  vector<int> array_strat_mean (500, 0);
  for (int i = 0; i != nind*ncomp; i++) {
    if (strat_mean[i] > 499) {cerr<<"No more than 500 strata are allowed\n";exit(1);}
    array_strat_mean[ strat_mean[i] ]++;
    if ( array_strat_mean[ strat_mean[i] ] == 1 ) nstrat_mean++;
  }

  int nstrat_assoc = 1;
  //vector<int> array_strat_assoc (500, 0);
  //for (int i = 0; i != nind*ncomp; i++) {
  //  if (strat_assoc[i] > 499) {cerr<<"No more than 500 strata are allowed\n";exit(1);}
  //  array_strat_assoc[ strat_assoc[i] ]++;
  //  if ( array_strat_assoc[ strat_assoc[i] ] == 1 ) nstrat_assoc++;
  //}
  //nstrat_assoc = 1;  ///here I override the sue of stratification in the logistic regression

  //if (nstrat_var > 1) {cout<<"Using stratification for variances: "<< nstrat_var <<" strata\n";}


  // Which frequency model
  MODEL fit_model = DISEASE; 
  if(  *INTEGER(pi_model) == 0 ) fit_model = DISEASE;
  else if(  *INTEGER(pi_model) == 1 ||  length(design_disease) <= 1 ) fit_model = HETERO;
  else if( *INTEGER(pi_model) == 2 ) fit_model = QT;
  else{
     Rprintf("Invalid model parameter : %d \n",*INTEGER(pi_model) );
  }

  // Which hypothesis
  HYPOTHESIS h = (test == "H0") ? H0 : H1;

  // Component type
  int component_type = *INTEGER(mix_model) <= 20 ? 1 : 2;

  // Constraint type
  //int signal_model = *INTEGER(mix_model) - component_type*10;

  
  // Dont update means/variances if there is no design matrix
  //bool fix_means = length(design_mean) > 1 ? false : true;
  //bool fix_vars = length(design_variance) > 1 ? false : true;
  CNV_signal * myCNV = new CNV_signal(nind, ncomp, cohort, signal, disease, 
				      mean_start, var_start, nu_start, alpha_start, 
				      logit_offset_p,
				      mean_design, var_design, disease_design, 
				      mean_dims[1], variance_dims[1], disease_dims[1], 
				      fit_model, h, logP_threshold, min_freq, 
				      strat_var, nstrat_var, strat_mean, nstrat_mean, strat_assoc, nstrat_assoc);
  
  string status;
  vector<double> dpost;


  // Fit the required model
  if(component_type == 1){
    int signal_model = *INTEGER(mix_model) - 10;
    fit_model_gaussian( myCNV, dpost, status, max_iter, tol, signal_model);
  }
  if(component_type == 2){
    int signal_model = *INTEGER(mix_model) - 2000;
    fit_model_t( myCNV, dpost, status, max_iter, tol, signal_model );
  }

  //Interface with R
  SEXP post, stat, ret;

  PROTECT(post = allocMatrix(REALSXP,nind*ncomp, 8));
  double * dres = REAL(post);
  //vector<double> dpost = myCNV->GetPosterior();       // Get the final data frame 
  for (size_t i = 0; i != dpost.size(); i++)   {
    dres[i] = dpost[i]; 
  }

  PROTECT(stat=allocVector(STRSXP,1));                // Allocate storage for status string
  SET_STRING_ELT(stat, 0, mkChar( status.c_str() ));

  PROTECT( ret = allocVector(VECSXP, 2) );            // Allocate and fill return list 
  SET_VECTOR_ELT(ret, 0, post);
  SET_VECTOR_ELT(ret, 1, stat);

  delete myCNV;

  UNPROTECT(3);
  return ret;
}

void fit_model_t(CNV_signal* model, vector<double>& posterior, string& status, const double& max_iter, const double& tol, int& model_type)
{
  //Work out the model
  int mean_flag = (model_type - model_type%100)/100;
  model_type -= mean_flag*100;
  int var_flag = (model_type - model_type%10)/10;
  model_type -= var_flag*10;
  int nu_flag = model_type;

  if( mean_flag != 3 && mean_flag != 4) {
    cerr << "error in mean_flag : " << mean_flag << ". aborting." << endl; exit(1);
  }
  if( var_flag != 1 && var_flag != 2 && var_flag != 3 && var_flag != 4 ){
    cerr << "error in var_flag : " << var_flag << ". aborting." << endl; exit(1);
  }
  if( nu_flag != 1 && nu_flag != 2 && nu_flag != 3 && nu_flag != 4 ){
    cerr << "error in nu_flag : " << nu_flag << ". aborting." << endl; exit(1);
  }
  //cout << mean_flag << "\t" << var_flag << "\t" << nu_flag << endl;

  double oldLike = 0;
  vector<double> old_posterior;
  try{
    model->ExpectationT();            // Expectation step
    model->ComputePosterior();                       // Computes the posterior distribution
    oldLike = model->GetLogLikelihood();
    old_posterior = model->GetPosterior();
  }
  catch(MyException e){
    status = "F"; 
    cerr << e.message << endl;
  }
  
  int loop = 0;
  while (1) {

    double newLike = 0;
    try{
            
      model->ExpectationT();                                      
      model->ComputePosterior();    
      
      model->MaximizeAlpha();
      model->MaximizeMeansT(mean_flag);  
      model->MaximizeVariancesT(var_flag);   
      model->MaximizeNuT(nu_flag);
      model->Check_order();
      
      model->ExpectationT(); 
      model->ComputePosterior();

      newLike = model->GetLogLikelihood();
    }
    catch(MyException e){
      status = "F"; 
      cerr << e.message << endl;
      break;
    }
 
    if (loop % 10 == 0) {       //////// final convergence test, every 10 iterations
      //model->PrintParams();
      //cout << "lnL : " << newLike << endl;
      if (fabs(newLike - oldLike) < tol){
	status = "C";
	break;
      }      
      swap(newLike, oldLike);
      
    }
    if (loop == max_iter){
      //Do a final check since we may have reached convergence
      if (fabs(newLike - oldLike) < tol){
	status = "C";
	break;
      } 
      else{
	status = "M";
	break;
      }    
    }  

    loop++;
  }
  //cout << "niter : " << loop << "\t" << status << endl;

  posterior = model->GetPosterior(); 
}
  
void fit_model_gaussian(CNV_signal* model, vector<double>& posterior, string& status,const double& max_iter, const double& tol, int& model_type)
{ 
  double oldLike = 0;
  try{
    model->ExpectationG();            // Expectation step
    model->ComputePosterior();                       // Computes the posterior distribution
    oldLike = model->GetLogLikelihood();
  }
  catch(MyException e){
    status = "F"; 
    cerr << e.message << endl;
  }
  
  int loop = 0;
  while (1) {
    //cout<<"Iter "<<loop<<endl;
    double newLike = 0;
    try{
      
      model->MaximizeMeansG();                                  // Maximization step for means 

      if(model_type == 0) model->MaximizeVariancesG();               // Maximization step for variances 
      else model->MaximizeVariancesPosteriorG(model_type);           // Maximization step for variances using posterior distribution L(data)*prior    
      

      model->Check_order();
      model->ExpectationG();                                      // Expectation step
      model->ComputePosterior();                                  // Computes the posterior distribution
      model->MaximizeAlpha();                                     // Maximization for alpha, disease, qt
      if (model->proba_disease[0] != model->proba_disease[0]) exit(1);

      model->ExpectationG();                                      // Expectation step
      model->ComputePosterior();                                  // Computes the posterior distribution
      newLike = model->GetLogLikelihood();
      //cout<<"aa "<<model->mean[0]<<"\t"<<model->mean[1100]<<"\t"<<model->mean[2200]<<endl;

    }
    catch(MyException e){
      status = "F"; 
      cerr << e.message << endl;
      break;
    }
    
    if (loop % 10 == 0) {       //////// final convergence test, every 10 iterations
      //cout << "lnL : " << newLike << endl;
      if (fabs(newLike - oldLike) < tol){
	status = "C";
	break;
      }      
      swap(newLike, oldLike);
      
    }
    loop++;
    if (loop == max_iter){
      //Do a final check since we may have reached convergence
      if (fabs(newLike - oldLike) < tol){
	status = "C";
	break;
      } 
      else{
	status = "M";
	break;
      }    
    }  
  }

  posterior = model->GetPosterior(); 
}



SEXP getListElement(SEXP list, const char *str)    
{      
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);    
  int i;
  for (i = 0; i < length(list); i++)        
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {      
      elmt = VECTOR_ELT(list, i);           
      break;        
    }      
  return elmt;     
}  



