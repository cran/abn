#include <R.h>
#include <Rdefines.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "structs.h"
#include "utility.h" 
#include "network_score.h"
#include <gsl/gsl_matrix.h>

#define DEBUG_12

                                         
SEXP fitabn_marginals(SEXP R_obsdata, SEXP R_dag,SEXP R_numVars, SEXP R_vartype, SEXP R_maxparents,SEXP R_priors_mean, SEXP R_priors_sd,SEXP R_priors_gamshape,SEXP R_priors_gamscale,
		      SEXP R_maxiters, SEXP R_epsabs, SEXP R_verbose, SEXP R_errorverbose, SEXP R_groupedvars, SEXP R_groupids, SEXP R_epsabs_inner,SEXP R_maxiters_inner,
	              SEXP R_finitestepsize, SEXP R_hparams,
		      SEXP R_childid, SEXP R_paramid, SEXP R_denom_modes, SEXP R_betafixed, SEXP R_mlik, SEXP R_maxiters_hessian)
{
/** ****************/
/** declarations **/

int errverbose,verbose;
datamatrix data,designmatrix;
const double priormean=asReal(R_priors_mean);/*Rprintf("priormean=%f %f\n",priormean[0],priormean[5]);*/
const double priorsd=asReal(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double priorgamshape=asReal(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double priorgamscale=asReal(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
/*const int *vartype=INTEGER(R_vartype);*/
/*const int numvariates=asInteger(R_numvariates);*/
int numVars=asInteger(R_numVars);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%u ",vartype[i]);} Rprintf("\n");*/
network dag;
/*cycle cyclestore;*/
/*storage nodescore;*/
SEXP posterior;
double *denom_modes=REAL(R_denom_modes);
int childid=asInteger(R_childid);
int paramid=asInteger(R_paramid);

/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
const int maxiters=asInteger(R_maxiters);
const double epsabs=asReal(R_epsabs);
int maxiters_inner=asInteger(R_maxiters_inner);
int maxiters_hessian=asInteger(R_maxiters_hessian);
double epsabs_inner=asReal(R_epsabs_inner);
const int maxparents=asInteger(R_maxparents);
double finitestepsize=asReal(R_finitestepsize);
/*double h_lowerend=REAL(R_hparams)[0];
double h_upperend=REAL(R_hparams)[1];*/
double h_guess=REAL(R_hparams)[0];
double h_epsabs=REAL(R_hparams)[1];
double betafixed=asReal(R_betafixed);
double mlik=asReal(R_mlik);
verbose=asInteger(R_verbose);
errverbose=asInteger(R_errorverbose);
/** end of argument parsing **/

/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
#ifdef JUNK
PROTECT(listresults = allocVector(VECSXP, 2));
for(i=0;i<2;i++){
				PROTECT(tmplistentry=NEW_NUMERIC(numvariates));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
#endif	

							
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */



make_dag(&dag, numVars,R_dag,0,R_vartype,&maxparents,R_groupedvars);/** user supplied DAG **/
  /*printDAG(&dag,2);*/
make_data(R_obsdata,&data,R_groupids);
#ifdef JUNK
/** An R MATRIX it is single dimension and just needs to be unrolled */
posterior=gsl_matrix_alloc(numvariates,2);
for(j=0;j<2;j++){for(i=0;i<numvariates;i++){gsl_matrix_set(posterior,i,j,REAL(R_posterior)[i+j*numvariates]);}} 
#endif
/*gsl_set_error_handler_off();*//*Rprintf("Warning: turning off GSL Error handler\n"); */

/*calc_network_Score(&dag,&data,&designmatrix,
		   priormean,priorsd,priorgamshape,priorgamscale,
		   maxiters,epsabs,verbose,errverbose,listresults,1,epsabs_inner,maxiters_inner,finitestepsize,
		   h_lowerend, h_upperend, h_guess, h_epsabs);*/
	  
PROTECT(posterior=NEW_NUMERIC(1));				
calc_parameter_marginal(&dag,&data,&designmatrix,
		   priormean,priorsd,priorgamshape,priorgamscale,
		   maxiters,epsabs,verbose,errverbose, 
		   denom_modes,childid,paramid,
			/*pdfminval, pdfstepsize,*/
			epsabs_inner,maxiters_inner,finitestepsize, h_guess, h_epsabs,maxiters_hessian,betafixed, mlik,REAL(posterior) );

/*gsl_set_error_handler (NULL);*//*Rprintf("Restoring: GSL Error handler\n");*/ /** restore the error handler*/

/** set up memory storage for any network which has each node with <=maxparents **/

/** Roll back into an R MATRIX */
/*for(j=0;j<2;j++){for(i=0;i<numvariates;i++){REAL(VECTOR_ELT(listresults,j))[i]=gsl_matrix_get(posterior,i,j);}} */

#ifdef JUNK
gsl_matrix_free(posterior);

UNPROTECT(1);

return(listresults);
#endif

 
UNPROTECT(1);
return(posterior);

}


