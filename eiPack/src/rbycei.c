

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "eiutil.h"



SEXP
rbycei_fcn1 (SEXP alphamatrix,
	   SEXP betaarray,
	   SEXP TT,
	   SEXP XX,
	   SEXP tuneA,
	   SEXP tuneB,
           SEXP NG,
	   SEXP NP,
	   SEXP Precincts,
	   SEXP Lambda1,
	   SEXP Lambda2,
	   SEXP Sample,
	   SEXP Thin,
	   SEXP Burnin,
	   SEXP Verbose,
	     SEXP Savebeta,
	     SEXP RR,
	     SEXP betanames
	     ){

  int nProtected = 0, nComps = 0, counter = 0;

  R_len_t ii, rr, cc, tt, qq, kk; 
  SEXP lbm,ccount, hldr; 
  SEXP a_acc, b_acc, a_draws, b_draws, ccount_draws, output_list;  
  R_len_t ng, np, prec, thin, samp, burn, iters, verbose;
  double lambda_1, lambda_2;
  double aprop, acurr, asumm1, lbm_rc, aprop_ll, acurr_ll;
  double bprop, bprop_ref, bcurr, bcurr_ref, tt_ci, tt_Ci, bprop_x, bprop_ref_x, bcurr_x, bcurr_ref_x, ulim, bcurr_ll, bprop_ll;

  ng = INTEGER(NG)[0];
  np = INTEGER(NP)[0];
  prec = INTEGER(Precincts)[0];
  lambda_1 = REAL(Lambda1)[0];
  lambda_2 = REAL(Lambda2)[0];

  samp = INTEGER(Sample)[0];
  thin = INTEGER(Thin)[0];
  burn = INTEGER(Burnin)[0];
  iters = burn + samp*thin;
  verbose = INTEGER(Verbose)[0];


  /*
   */

  PROTECT(hldr = allocVector(NILSXP, 1));
  ++nProtected;

  /*
  PROTECT(lbm = allocMatrix(REALSXP, ng , np));
  ++nProtected;
  */

  PROTECT(a_acc = allocVector(REALSXP, ng * np));
  ++nProtected;
  ++nComps;

  PROTECT(b_acc = allocVector(REALSXP, ng*(np-1)*prec));
  ++nProtected;
  ++nComps;


  PROTECT(ccount = allocMatrix(REALSXP, ng, np));
  ++nProtected;
  ++nComps;

  for(qq = 0; qq < ng*np; ++qq){
    REAL(a_acc)[qq] = 0; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = 0; }

  PROTECT(ccount_draws = allocMatrix(REALSXP, samp , ng*np));
  ++nProtected;
  ++nComps;


  PROTECT(a_draws = allocMatrix(REALSXP, samp , ng*np));
  ++nProtected;
  ++nComps;


  if(INTEGER(Savebeta)[0] == 0){
  PROTECT(b_draws = allocMatrix(REALSXP, samp, ng*np*prec));
  ++nProtected;
  ++nComps;
  }

  GetRNGstate();


  for(kk = 0; kk < iters; ++kk){
  for(ii = 0; ii < prec; ++ii){
    for(rr = 0; rr < ng; ++rr){
      for(cc = 0; cc < (np - 1); ++cc){
	bcurr = REAL(betaarray)[rr + ng*cc + ng*np*ii];
	bcurr_ref = REAL(betaarray)[rr + ng*(np-1) + ng*np*ii];
	ulim = bcurr + bcurr_ref;
	bprop = rnorm(bcurr, REAL(tuneB)[rr + ng*cc + ng*(np-1)*ii]);
	bprop_ref = ulim - bprop;


	if(bprop > 0 && bprop < ulim){
	  tt_ci = REAL(TT)[cc + np*ii];
	  tt_Ci = REAL(TT)[(np - 1) + np*ii];
	  
	  bcurr_x = 0; bcurr_ref_x = 0;
	  for(qq = 0; qq < ng; ++qq){
	    bcurr_x += REAL(betaarray)[qq + ng*cc + ng*np*ii]* REAL(XX)[qq + ng*ii];
	    bcurr_ref_x += REAL(betaarray)[qq + ng*(np-1) + ng*np*ii] * REAL(XX)[qq + ng*ii];
	  }
	 

	  bprop_x = bcurr_x - bcurr*REAL(XX)[rr + ng*ii] + bprop*REAL(XX)[rr + ng*ii];
	  bprop_ref_x =  bcurr_ref_x - bcurr_ref*REAL(XX)[rr + ng*ii] + bprop_ref*REAL(XX)[rr + ng*ii];

	  bprop_ll = beta_ll(bprop, bprop_ref, REAL(alphamatrix)[rr + ng*cc], REAL(alphamatrix)[rr + ng*(np-1)],tt_ci, tt_Ci, bprop_x, bprop_ref_x);
	  bcurr_ll = beta_ll(bcurr, bcurr_ref,  REAL(alphamatrix)[rr + ng*cc], REAL(alphamatrix)[rr + ng*(np-1)], tt_ci, tt_Ci, bcurr_x, bcurr_ref_x);
	
	   /* Rprintf("%f %f %f %f %f\n", bprop, bprop_ref, bprop_x, bprop_ref_x, bprop_ll - bcurr_ll); 
	    */
	  if(acc_tog(bprop_ll, bcurr_ll) == 1){
	    REAL(betaarray)[rr + ng*cc + ng*np*ii] = bprop;
	    REAL(betaarray)[rr + ng*(np-1) + ng*np*ii] = bprop_ref;
	    REAL(b_acc)[rr + ng*cc + ng*(np - 1)*ii] += 1;
  
	  }
	}
      }
    }
  }


  PROTECT(lbm = logbetamat(betaarray, NG, NP, Precincts));

  for(rr = 0; rr < ng; ++rr){
    for(cc = 0; cc < np; ++cc){
      
      acurr = REAL(alphamatrix)[rr + ng*cc];
      aprop = rnorm(acurr, REAL(tuneA)[rr + ng*cc]);
     
      lbm_rc = REAL(lbm)[rr + ng*cc];
      asumm1 = 0; 
      for(tt = 0; tt < np; ++tt){
	asumm1 += REAL(alphamatrix)[rr + ng*tt];
      }
      asumm1 = asumm1 - acurr;

      if(aprop > 0){
	aprop_ll = alpha_ll(aprop, asumm1, lbm_rc, lambda_1, lambda_2, prec);
	acurr_ll = alpha_ll(acurr, asumm1, lbm_rc, lambda_1, lambda_2, prec);
	/*Rprintf("%f %f \n", aprop, aprop_ll - acurr_ll);
	 */
	if(acc_tog(aprop_ll, acurr_ll) == 1){
	  REAL(alphamatrix)[rr + ng*cc] = aprop;
	  REAL(a_acc)[rr + ng*cc] += 1;
	 
	}
      }
    }
  }
  UNPROTECT(1);



  if(kk >= burn && ((kk % thin) == 0)){


    ccount = cellcount(betaarray,RR, NG, NP, Precincts);


    for(qq = 0; qq < np*ng; ++qq){
      REAL(a_draws)[counter + qq*samp] = REAL(alphamatrix)[qq];
      REAL(ccount_draws)[counter + qq*samp] = REAL(ccount)[qq];
	}

    if(INTEGER(Savebeta)[0] == 0){
     for(qq = 0; qq < np*ng*prec; ++qq){
       REAL(b_draws)[counter + qq*samp] = REAL(betaarray)[qq];
	}
    }


    if(INTEGER(Savebeta)[0] == 2){
      write_beta(betaarray, betanames);
    }

     counter += 1;
  }

  if(verbose > 0 && kk % verbose == 0){
    Rprintf("\n MCMC iteration %i of %i \n", kk + 1, iters); 
  }

  R_CheckUserInterrupt();

  }

   /*
    *   PROTECT(dim_matrix1 = allocVector(INTSXP, 2));
    *++nProtected; 
    *   INTEGER(dim_matrix1)[0] = ncol;
    *   INTEGER(dim_matrix1)[1] = nrow;
    *   
    *   setAttrib(matrix1, R_DimSymbol, dim_matrix1);
    */


  for(qq = 0; qq < ng*np; ++qq){
    REAL(a_acc)[qq] = REAL(a_acc)[qq]/iters; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = REAL(b_acc)[qq]/iters; }


  if(INTEGER(Savebeta)[0]==0){
    PROTECT(output_list = allocVector(VECSXP, 5));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, a_draws);
    SET_VECTOR_ELT(output_list, 1, b_draws);
    SET_VECTOR_ELT(output_list, 2, a_acc);    
    SET_VECTOR_ELT(output_list, 3, b_acc);
    SET_VECTOR_ELT(output_list, 4, ccount_draws);
  }else{
    PROTECT(output_list = allocVector(VECSXP, 4));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, a_draws);
    SET_VECTOR_ELT(output_list, 1, a_acc);    
    SET_VECTOR_ELT(output_list, 2, b_acc);
        SET_VECTOR_ELT(output_list, 3, ccount_draws);
  }


  PutRNGstate();
  UNPROTECT(nProtected); 

  return(output_list); 

}
