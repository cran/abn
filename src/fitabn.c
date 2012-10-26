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

#define DEBUG_12

                                         
SEXP fitabn(SEXP R_obsdata, SEXP R_dag,SEXP R_numVars, SEXP R_vartype, SEXP R_maxparents,SEXP R_priors_mean, SEXP R_priors_sd,SEXP R_priors_gamshape,SEXP R_priors_gamscale,
			 SEXP R_maxiters, SEXP R_epsabs, SEXP R_verbose, SEXP R_errorverbose, SEXP R_groupedvars, SEXP R_groupids, SEXP R_epsabs_inner,SEXP R_maxiters_inner,
	                 SEXP R_finitestepsize, SEXP R_hparams, SEXP R_maxitershessian)
{
/** ****************/
/** declarations **/
unsigned int i;
int errverbose,verbose;
datamatrix data,designmatrix;
const double priormean=asReal(R_priors_mean);/*Rprintf("priormean=%f %f\n",priormean[0],priormean[5]);*/
const double priorsd=asReal(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double priorgamshape=asReal(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double priorgamscale=asReal(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
/*const int *vartype=INTEGER(R_vartype);*/
int numVars=asInteger(R_numVars);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%u ",vartype[i]);} Rprintf("\n");*/
network dag;
/*cycle cyclestore;*/
/*storage nodescore;*/
SEXP listresults;
SEXP tmplistentry;
/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
const int maxiters=asInteger(R_maxiters);
const double epsabs=asReal(R_epsabs);
int maxiters_inner=asInteger(R_maxiters_inner);
int maxiters_hessian=asInteger(R_maxitershessian);
double epsabs_inner=asReal(R_epsabs_inner);
const int maxparents=asInteger(R_maxparents);
double finitestepsize=asReal(R_finitestepsize);
/*double h_lowerend=REAL(R_hparams)[0];
double h_upperend=REAL(R_hparams)[1];*/
double h_guess=REAL(R_hparams)[0];
double h_epsabs=REAL(R_hparams)[1];

verbose=asInteger(R_verbose);
errverbose=asInteger(R_errorverbose);
/*Rprintf("h params=%e %e %e %e\n",h_lowerend,h_upperend,h_guess,h_epsabs);*/
/** end of argument parsing **/

/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
PROTECT(listresults = allocVector(VECSXP, numVars));
for(i=0;i<numVars;i++){
				PROTECT(tmplistentry=NEW_NUMERIC(numVars+2+4));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */
/** Check for grouped variables e.g. random effects  */

  make_dag(&dag, numVars,R_dag,0,R_vartype,&maxparents,R_groupedvars);/** user supplied DAG **/
  /*printDAG(&dag,2);*/
  make_data(R_obsdata,&data,R_groupids);
  
  
/** some print debugging about grouping part **/
/*if(LENGTH(R_groupids)==1){Rprintf("no grouping %d %d\n", groupedvars[0],groupids[0]);
} else {Rprintf("variables to be grouped:\n");for(i=0;i<LENGTH(R_groupedvars);i++){Rprintf("%d",groupedvars[0]);}Rprintf("\n");
        Rprintf("Grouping IDs:\n");for(i=0;i<LENGTH(R_groupids);i++){Rprintf("%d ",groupids[i]);}Rprintf("\n");
}*/
/** *****/

/*gsl_set_error_handler_off();*//*Rprintf("Note: turning off GSL Error handler\n");*/ 

calc_network_Score(&dag,&data,&designmatrix,
		   priormean,priorsd,priorgamshape,priorgamscale,
		   maxiters,epsabs,verbose,errverbose,listresults,1,epsabs_inner,maxiters_inner,finitestepsize,
		    h_guess, h_epsabs, maxiters_hessian);
		   /** NOTE: last arg to calc_network_Score is to store modes - this does no harm
		   and is potentially useful is fitabn.R also wants posterior dists*/

/*gsl_set_error_handler (NULL);*//** restore the error handler*/

        /** set up memory storage for any network which has each node with <=maxparents **/

free_dag(&dag);/** free any memory not allocated using Ralloc() */
UNPROTECT(1);

return(listresults);

}


