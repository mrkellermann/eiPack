

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
rbycei_fcn4 (SEXP drvector,
	     SEXP betaarray,
	     SEXP gammamatrix,
	     SEXP deltamatrix,
	     SEXP TT,
	     SEXP XX,
	     SEXP ZZ,
	     SEXP tuneDr,
	     SEXP tuneB,
	     SEXP tuneG,
	     SEXP tuneD,
	     SEXP NG,
	     SEXP NP,
	     SEXP Precincts,
	     SEXP Lambda1,
	     SEXP Lambda2,
	     SEXP Covprior, 
	     SEXP Delmean,
	     SEXP Delsd,
	     SEXP Gammean,
	     SEXP Gamsd,
	     SEXP Sample,
	     SEXP Thin,
	     SEXP Burnin,
	     SEXP Verbose,
	     SEXP Savebeta,
	     SEXP RR,
	     SEXP usr_fcn, 
	     SEXP usr_env,
	     SEXP usr_len,
	     SEXP betanames
	     ){

  int nProtected = 0, nComps = 0, counter = 0;

  R_len_t ii, rr, cc, qq, kk; 
  SEXP ccount, usr, hldr; 
  SEXP usr_list, dr_dim, beta_dim, gamma_dim, delta_dim, TT_dim, RR_dim; 
  SEXP dr_acc, b_acc, g_acc, d_acc; 
  SEXP dr_draws, b_draws, g_draws, d_draws, ccount_draws, usr_draws, output_list;  
  R_len_t ng, np, npm1, prec, thin, samp, burn, iters, verbose,usr_length;
  double lambda1, lambda2;
  SEXP explp;
  double bprop, bprop_ref, bcurr, bcurr_ref, tt_ci, tt_Ci, bprop_x, bprop_ref_x, bcurr_x, bcurr_ref_x, ulim, bcurr_ll, bprop_ll, alpha_ci, alpha_Ci;
  double drcurr, drprop, drcurr_ll, drprop_ll, tmp_currB, tmp_propB;
  double gcurr, gprop, gcurr_ll, gprop_ll, tmp_gcurr, tmp_gprop;
  double dcurr, dprop, dcurr_ll, dprop_ll, tmp_dcurr, tmp_dprop;
  int covprior;

  ng = INTEGER(NG)[0];
  np = INTEGER(NP)[0];
  prec = INTEGER(Precincts)[0];
  npm1 = np - 1;
  lambda1 = REAL(Lambda1)[0];
  lambda2 = REAL(Lambda2)[0];
  covprior = INTEGER(Covprior)[0];  

  samp = INTEGER(Sample)[0];
  thin = INTEGER(Thin)[0];
  burn = INTEGER(Burnin)[0];
  iters = burn + samp*thin;
  verbose = INTEGER(Verbose)[0];
  usr_length = INTEGER(usr_len)[0];



  PROTECT(beta_dim = allocVector(INTSXP, 3));
  ++nProtected; 
  INTEGER(beta_dim)[0] = ng;
  INTEGER(beta_dim)[1] = np;
  INTEGER(beta_dim)[2] = prec;
  setAttrib(betaarray, R_DimSymbol, beta_dim);

  PROTECT(gamma_dim = allocVector(INTSXP, 2));
  ++nProtected; 
  INTEGER(gamma_dim)[0] = ng;
  INTEGER(gamma_dim)[1] = np;
  setAttrib(gammamatrix, R_DimSymbol, gamma_dim);

  PROTECT(delta_dim = allocVector(INTSXP, 2));
  ++nProtected; 
  INTEGER(delta_dim)[0] = ng;
  INTEGER(delta_dim)[1] = np;
  setAttrib(deltamatrix, R_DimSymbol, delta_dim);

  PROTECT(TT_dim = allocVector(INTSXP, 2));
  ++nProtected; 
  INTEGER(TT_dim)[0] = np;
  INTEGER(TT_dim)[1] = prec;
  setAttrib(TT, R_DimSymbol, TT_dim);

  PROTECT(RR_dim = allocVector(INTSXP, 2));
  ++nProtected; 
  INTEGER(RR_dim)[0] = ng;
  INTEGER(RR_dim)[1] = prec;
  setAttrib(RR, R_DimSymbol, RR_dim);

  


  /*
   */
  
  PROTECT(hldr = allocVector(NILSXP, 1));
  ++nProtected;
  PROTECT(ccount = allocMatrix(REALSXP, ng, np));
  ++nProtected;
  PROTECT(explp = allocVector(REALSXP, ng*np*prec));
  ++nProtected;

  PROTECT(dr_acc = allocVector(REALSXP, ng));
  ++nProtected;
  ++nComps;
  PROTECT(b_acc = allocVector(REALSXP, ng*(np-1)*prec));
  ++nProtected;
  ++nComps;
  PROTECT(g_acc = allocVector(REALSXP, ng*(np-1)));
  ++nProtected;
  ++nComps;
  PROTECT(d_acc = allocVector(REALSXP, ng*(np-1)));
  ++nProtected;
  ++nComps;



  for(qq = 0; qq < ng; ++qq){
    REAL(dr_acc)[qq] = 0; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = 0; }
  for(qq = 0; qq < ng*(np-1);++qq){
    REAL(g_acc)[qq] = 0;
    REAL(d_acc)[qq] = 0;}

  PROTECT(usr = allocVector(REALSXP, usr_length));
  ++nProtected;
  ++nComps;
  PROTECT(ccount_draws = allocMatrix(REALSXP, samp , ng*np));
  ++nProtected;
  ++nComps;

  PROTECT(usr_list = allocVector(VECSXP, 6));
  ++nProtected;

  PROTECT(usr_draws = allocMatrix(REALSXP, samp , usr_length));
  ++nProtected;
  ++nComps;

  PROTECT(dr_draws = allocMatrix(REALSXP, samp , ng));
  ++nProtected;
  ++nComps;

  PROTECT(g_draws = allocMatrix(REALSXP, samp , ng*npm1));
  ++nProtected;
  ++nComps;

  PROTECT(d_draws = allocMatrix(REALSXP, samp , ng*npm1));
  ++nProtected;
  ++nComps;

  if(INTEGER(Savebeta)[0] == 0){
  PROTECT(b_draws = allocMatrix(REALSXP, samp, ng*np*prec));
  ++nProtected;
  ++nComps;
  }


  GetRNGstate();

  for(kk = 0; kk < iters; ++kk){

    /* create exp(gamma_rc + delta_rc*Z_i)*/

    for(ii = 0; ii < prec; ++ii){
      for(rr = 0; rr < ng; ++rr){
	for(cc = 0; cc < np; ++cc){
	  REAL(explp)[rr + ng*cc + ng*np*ii] = exp(REAL(gammamatrix)[rr + ng*cc] + REAL(deltamatrix)[rr + ng*cc]*REAL(ZZ)[ii]);
	}
      }
    }


    /* draw d_r */

    for(rr = 0; rr < ng; ++rr){
      drcurr = REAL(drvector)[rr];
      drprop = rnorm(drcurr, REAL(tuneDr)[rr]);
      if(drprop > 0){
	drcurr_ll = 0;
	drprop_ll = 0;
	for(ii = 0; ii < prec; ++ii){
	  for(cc = 0; cc < np; ++cc){
	    drcurr_ll += -1*lgammafn(drcurr * REAL(explp)[rr + ng*cc + ng*np*ii]) + drcurr*REAL(explp)[rr + ng*cc + ng*np*ii]*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
	    drprop_ll += -1*lgammafn(drprop * REAL(explp)[rr + ng*cc + ng*np*ii]) + drprop*REAL(explp)[rr + ng*cc + ng*np*ii]*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
	  }
	}
	for(ii = 0; ii < prec; ++ii){
	  tmp_currB = 0;
	  tmp_propB = 0;
	  for(cc = 0; cc < np; ++cc){
	    tmp_currB += drcurr * REAL(explp)[rr + ng*cc + ng*np*ii];
	    tmp_propB += drprop * REAL(explp)[rr + ng*cc + ng*np*ii];
	  }
	  drcurr_ll += lgammafn(tmp_currB);
	  drprop_ll += lgammafn(tmp_propB);
	}
      
	drcurr_ll = drcurr_ll - lambda2*drcurr + (lambda1 - 1)*log(drcurr);
	drprop_ll = drprop_ll - lambda2*drprop + (lambda1 - 1)*log(drprop);


	/*	Rprintf("%f %f %f %f\n", drcurr, drprop, drprop_ll, drcurr_ll); */

      if(acc_tog(drprop_ll, drcurr_ll) == 1){
	REAL(drvector)[rr] = drprop;
	REAL(dr_acc)[rr] += 1;
      }
      }
    }

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

	  alpha_ci = REAL(drvector)[rr]*REAL(explp)[rr + ng*cc + ng*np*ii];
	  alpha_Ci = REAL(drvector)[rr]*REAL(explp)[rr + ng*npm1 + ng*np*ii];

	  bprop_ll = beta_ll(bprop, bprop_ref, alpha_ci, alpha_Ci,tt_ci, tt_Ci, bprop_x, bprop_ref_x);
	  bcurr_ll = beta_ll(bcurr, bcurr_ref,  alpha_ci, alpha_Ci, tt_ci, tt_Ci, bcurr_x, bcurr_ref_x);
	
	   /*  
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



    

  for(rr = 0; rr < ng; ++rr){
    for(cc = 0; cc < npm1; ++cc){
      gcurr = REAL(gammamatrix)[rr + ng*cc];
      gprop = rnorm(gcurr, REAL(tuneG)[rr + cc*ng]);
      
      dcurr = REAL(deltamatrix)[rr + ng*cc];
      dprop = rnorm(dcurr, REAL(tuneD)[rr + cc*ng]);
      
      gcurr_ll = 0;
      gprop_ll = 0;
      for(ii = 0; ii < prec; ++ii){
	tmp_gcurr = exp(gcurr + dcurr*REAL(ZZ)[ii]);
	tmp_gprop = exp(gprop + dprop*REAL(ZZ)[ii]);
	gcurr_ll += -1*lgammafn(REAL(drvector)[rr] * tmp_gcurr) + REAL(drvector)[rr]*tmp_gcurr*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
	gprop_ll += -1*lgammafn(REAL(drvector)[rr] * tmp_gprop) + REAL(drvector)[rr]*tmp_gprop*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
      }
      
      
      /*	Rprintf("%f %f %f %f\n", gcurr, gprop, gprop_ll, gcurr_ll); */
      
      for(ii = 0; ii < prec; ++ii){
	tmp_currB = 0;
	tmp_propB = 0;
	for(qq = 0; qq < np; ++qq){
	  tmp_currB += REAL(drvector)[rr] * exp(REAL(gammamatrix)[rr + ng*qq] + REAL(deltamatrix)[rr + ng*qq]*REAL(ZZ)[ii]);
	  /* Rprintf("%i %f \n", qq,  exp(REAL(gammamatrix)[rr + ng*qq] + REAL(deltamatrix)[rr + ng*qq]*REAL(ZZ)[ii]));*/
	}
	tmp_propB = tmp_currB - REAL(drvector)[rr] * exp(gcurr + dcurr*REAL(ZZ)[ii]) + REAL(drvector)[rr] * exp(gprop + dprop*REAL(ZZ)[ii]);
	/*Rprintf("%i %f \n",  cc, exp(gcurr + REAL(deltamatrix)[rr + ng*cc]*REAL(ZZ)[ii]));*/
	gcurr_ll += lgammafn(tmp_currB);
	gprop_ll += lgammafn(tmp_propB);
	
	if(verbose > 0 && kk % verbose == 0){
	  /*    Rprintf("%f %f %f %f %f %f \n", gcurr, gprop, gcurr_ll, gprop_ll, tmp_currB, tmp_propB); */
	}
	
      }
      
      if(covprior==1){
	gcurr_ll += dnorm(gcurr,REAL(Gammean)[rr + ng*cc] ,REAL(Gamsd)[rr + ng*cc] , 1) + dnorm(dcurr,REAL(Delmean)[rr + ng*cc] , REAL(Delsd)[rr + ng*cc], 1);
	gprop_ll += dnorm(gprop,REAL(Gammean)[rr + ng*cc] ,REAL(Gamsd)[rr + ng*cc] , 1) + dnorm(dprop,REAL(Delmean)[rr + ng*cc] ,REAL(Delsd)[rr + ng*cc] , 1);
      }
      
      if(acc_tog(gprop_ll, gcurr_ll) == 1){
	REAL(gammamatrix)[rr + ng*cc] = gprop;
	REAL(deltamatrix)[rr + ng*cc] = dprop;
	REAL(g_acc)[rr + ng*cc] += 1;
      }
    }
  }
  

    /*    for(rr = 0; rr < ng; ++rr){
     *   for(cc = 0; cc < npm1; ++cc){
     *  dcurr = REAL(deltamatrix)[rr + ng*cc];
     *	dprop = rnorm(dcurr, REAL(tuneD)[rr + cc*ng]);
     *
     *	dcurr_ll = 0;
     *	dprop_ll = 0;
     *	for(ii = 0; ii < prec; ++ii){
     *	  tmp_dcurr = exp(REAL(gammamatrix)[rr + ng*cc] + dcurr*REAL(ZZ)[ii]);
     *	  tmp_dprop = exp(REAL(gammamatrix)[rr + ng*cc] + dprop*REAL(ZZ)[ii]);
     *	  dcurr_ll += -1*lgammafn(REAL(drvector)[rr] * tmp_dcurr) + REAL(drvector)[rr]*tmp_dcurr*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
     *	  dprop_ll += -1*lgammafn(REAL(drvector)[rr] * tmp_dprop) + REAL(drvector)[rr]*tmp_dprop*log(REAL(betaarray)[rr + ng*cc + ng*np*ii]);
     *	}
     *	for(ii = 0; ii < prec; ++ii){
     *	  tmp_currB = 0;
     *	  tmp_propB = 0;
     *	  for(qq = 0; qq < np; ++qq){
     *	    tmp_currB += REAL(drvector)[rr] * exp(REAL(gammamatrix)[rr + ng*qq] + REAL(deltamatrix)[rr + ng*qq]*REAL(ZZ)[ii]);
     *	  }
     *	  tmp_propB = tmp_currB - REAL(drvector)[rr] * exp(REAL(gammamatrix)[rr + ng*cc] + dcurr*REAL(ZZ)[ii]) + REAL(drvector)[rr] * exp(REAL(gammamatrix)[rr + ng*cc] + dprop*REAL(ZZ)[ii]);
     *	  dcurr_ll += lgammafn(tmp_currB);
     *	  dprop_ll += lgammafn(tmp_propB);
     *	}
     * 
     *	
     * if(acc_tog(dprop_ll, dcurr_ll) == 1){
     *	REAL(deltamatrix)[rr + ng*cc] = dprop;
     *	REAL(d_acc)[rr + ng*cc] += 1;
     * }
     * }
     *}
     */
  

  if(kk >= burn && ((kk % thin) == 0)){
  
    SET_VECTOR_ELT(usr_list, 0, drvector);
    SET_VECTOR_ELT(usr_list, 1, betaarray);
    SET_VECTOR_ELT(usr_list, 2, gammamatrix);
    SET_VECTOR_ELT(usr_list, 3, deltamatrix);
    SET_VECTOR_ELT(usr_list, 4, TT);
    SET_VECTOR_ELT(usr_list, 5, RR);

    ccount = cellcount(betaarray, RR, NG, NP, Precincts);

    usr = usr_fun_eval(usr_fcn, usr_list, usr_env, usr_length);


    for(qq = 0; qq < usr_length; ++qq){
      REAL(usr_draws)[counter + qq*samp] = REAL(usr)[qq];
	}

    for(qq = 0; qq < np*ng; ++qq){
      REAL(ccount_draws)[counter + qq*samp] = REAL(ccount)[qq];
	}

    for(qq = 0; qq < ng; ++qq){
      REAL(dr_draws)[counter + qq*samp] = REAL(drvector)[qq];
    }

    for(qq = 0; qq < ng*npm1; ++qq){
      REAL(g_draws)[counter + qq*samp] = REAL(gammamatrix)[qq];
      REAL(d_draws)[counter + qq*samp] = REAL(deltamatrix)[qq];
    }

    if(INTEGER(Savebeta)[0] == 0){
     for(qq = 0; qq < np*ng*prec; ++qq){
       REAL(b_draws)[counter + qq*samp] = REAL(betaarray)[qq];
	}
    }


    if(INTEGER(Savebeta)[0] == 2){
      hldr = write_beta(betaarray, betanames);
    }


     counter += 1;
  }

  if(verbose > 0 && kk % verbose == 0){
    Rprintf("\n MCMC iteration %i of %i \n", kk + 1, iters); 
  }

  R_CheckUserInterrupt();

  }

 


  for(qq = 0; qq < ng; ++qq){
    REAL(dr_acc)[qq] = REAL(dr_acc)[qq]/iters; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = REAL(b_acc)[qq]/iters; }
  for(qq = 0; qq < ng*(np-1);++qq){
    REAL(g_acc)[qq] = REAL(g_acc)[qq]/iters; } 
  for(qq = 0; qq < ng*(np-1);++qq){
    REAL(d_acc)[qq] = REAL(d_acc)[qq]/iters; }

  if(INTEGER(Savebeta)[0]==0){
    PROTECT(output_list = allocVector(VECSXP, 10));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, dr_draws);
    SET_VECTOR_ELT(output_list, 1, b_draws);
    SET_VECTOR_ELT(output_list, 2, g_draws);
    SET_VECTOR_ELT(output_list, 3, d_draws);
    SET_VECTOR_ELT(output_list, 4, dr_acc);    
    SET_VECTOR_ELT(output_list, 5, b_acc);
    SET_VECTOR_ELT(output_list, 6, g_acc);
    SET_VECTOR_ELT(output_list, 7, d_acc);
    SET_VECTOR_ELT(output_list, 8, ccount_draws);
    SET_VECTOR_ELT(output_list, 9, usr_draws);
  }


  if(INTEGER(Savebeta)[0]!=0){
    PROTECT(output_list = allocVector(VECSXP, 9));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, dr_draws);
    SET_VECTOR_ELT(output_list, 1, g_draws);
    SET_VECTOR_ELT(output_list, 2, d_draws);
    SET_VECTOR_ELT(output_list, 3, dr_acc);    
    SET_VECTOR_ELT(output_list, 4, b_acc);
    SET_VECTOR_ELT(output_list, 5, g_acc);
    SET_VECTOR_ELT(output_list, 6, d_acc);
    SET_VECTOR_ELT(output_list, 7, ccount_draws);
    SET_VECTOR_ELT(output_list, 8, usr_draws);
  

  }


  PutRNGstate();
  UNPROTECT(nProtected); 

  return(output_list); 

}
