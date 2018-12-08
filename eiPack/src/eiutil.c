

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "eiutil.h"


double beta_ll (double beta, double beta_ref, double alpha, double alpha_ref, double tt_ci, double tt_Ci, double beta_X, double beta_ref_X){

  double ll;
  ll = tt_ci*log(beta_X) + tt_Ci*log(beta_ref_X) + (alpha-1)*log(beta) + (alpha_ref - 1)*log(beta_ref);

  return(ll);
}





double alpha_ll (double alpha, double sum_alphaminus, double lbmat, double lambda1, double lambda2, int prec){

  double ll;
  ll = (lambda1 - 1)*log(alpha) - lambda2*alpha + prec*(lgammafn((alpha + sum_alphaminus)) - lgammafn(alpha)) + lbmat*alpha;

  return(ll);
}





int acc_tog (double prop_ll, double curr_ll){

  int toggle = 0;
  double log_rat; 
  log_rat = prop_ll - curr_ll;
  if(log(runif(0,1)) < log_rat){
    toggle = 1;}
  return(toggle);
}




SEXP
cellcount (SEXP betas,
	   SEXP RR,
       SEXP NG,
       SEXP NP,
       SEXP Precincts){

  int nProtect = 0; 
  R_len_t rr, cc, ii;
  int ng, np, prec;
  double tmp = 0.0;  
  
  ng = INTEGER(NG)[0];
  np = INTEGER(NP)[0];
  prec = INTEGER(Precincts)[0];

  SEXP ret_val; 
  PROTECT(ret_val = allocMatrix(REALSXP, ng, np));
  ++nProtect;

  for(rr = 0; rr < ng; ++rr){
    for(cc = 0; cc < np; ++cc){
  tmp = 0.0;
  for(ii = 0; ii < prec; ++ii){
    tmp += REAL(betas)[rr + ng * cc + ng * np * ii]*REAL(RR)[rr + ng*ii];
  }
  REAL(ret_val)[rr + ng * cc] = tmp;
    }
  }

  UNPROTECT(1);
  return(ret_val);
}



SEXP
logbetamat (SEXP aa, 
       SEXP NG,
       SEXP NP,
       SEXP Precincts){

  int nProtect = 0; 
  R_len_t rr, cc, ii;
  int ng, np, prec;
  double tmp = 0.0;  
  
  ng = INTEGER(NG)[0];
  np = INTEGER(NP)[0];
  prec = INTEGER(Precincts)[0];

  SEXP ret_val; 
  PROTECT(ret_val = allocMatrix(REALSXP, ng, np));
  ++nProtect;

  for(rr = 0; rr < ng; ++rr){
    for(cc = 0; cc < np; ++cc){
  tmp = 0.0;
  for(ii = 0; ii < prec; ++ii){
    tmp += log(REAL(aa)[rr + ng * cc + ng * np * ii]);
  }
  REAL(ret_val)[rr + ng * cc] = tmp;
    }
  }

  UNPROTECT(1);
  return(ret_val);
}



SEXP
write_beta (SEXP betaarray,
	    SEXP filenames){

  R_len_t nn, ii; 
  nn = length(filenames);
  for(ii = 0; ii < nn; ++ii){
    char tmp[500];
    sprintf(tmp, "echo \"%.16f\" | gzip >>  %s", REAL(betaarray)[ii], CHAR(STRING_ELT(filenames,ii)));
    system(tmp);
  }


  R_CheckUserInterrupt();

  return(R_NilValue);
  }






SEXP
usr_fun_eval(SEXP usr_fun,
	      SEXP cur_values,
	      SEXP usr_env,
	      int usr_len){

  SEXP R_fcall, usr_fcn_output;
  if(!isFunction(usr_fun)) error("`usr_fun' must be a function");
  if(!isEnvironment(usr_env)) error("`usr_env' must be an environment");
  PROTECT(R_fcall = lang2(usr_fun, R_NilValue));
  SETCADR(R_fcall, cur_values);
  PROTECT(usr_fcn_output = allocVector(REALSXP, usr_len));
  usr_fcn_output = eval(R_fcall, usr_env);
  UNPROTECT(2);
  return(usr_fcn_output);
}


