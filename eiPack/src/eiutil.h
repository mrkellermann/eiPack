

#ifndef EIUTIL_H
#define EIUTIL_H

double beta_ll (double beta, double beta_ref, double alpha, double alpha_ref, double tt_ci, double tt_Ci, double beta_X, double beta_ref_X);

double alpha_ll (double alpha, double sum_alphaminus, double lbmat, double lambda1, double lambda2, int prec);

int acc_tog (double prop_ll, double curr_ll);

SEXP cellcount (SEXP betas,
	   SEXP RR,
       SEXP NG,
       SEXP NP,
	   SEXP Precincts);

SEXP logbetamat (SEXP aa, 
                 SEXP NG,
                 SEXP NP,
		 SEXP Precincts);


SEXP write_beta (SEXP betaarray,
		 SEXP filenames);



SEXP
usr_fun_eval(SEXP usr_fun,
	      SEXP cur_values,
	      SEXP usr_env,
	     int usr_len);


#endif


