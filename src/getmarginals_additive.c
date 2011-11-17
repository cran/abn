#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "utility_fns.h"
#include "network.h"
#include "network_laplace.h"

#include "laplace.h"
#include "laplace_marginals.h"
#include <gsl/gsl_errno.h>
#define DEBUG_12

SEXP getmarginals_additive(SEXP R_obsdata, SEXP R_dag,SEXP R_priors_mean, SEXP R_priors_sd, SEXP R_priors_gamshape,SEXP R_priors_gamscale,
			                     SEXP R_maxparents, SEXP R_verbose, SEXP R_vartype,
                           SEXP R_posterior, SEXP R_numvariates,SEXP R_whichnode, SEXP R_whichvariable, SEXP R_whichgaus, SEXP R_maxiters, SEXP R_epsabs)
{
/** ****************/
/** declarations **/
unsigned int i,j,maxparents;
unsigned int verbose;
unsigned int /*marginals,*/ numvariates,whichnode,whichvariable,whichgaus;
datamatrix obsdata, designmatrix;
const double *priormean=REAL(R_priors_mean);/*Rprintf("priormean=%f\n",priormean[0],priormean[5]);*/
const double *priorsd=REAL(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double *priorgamshape=REAL(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double *priorgamscale=REAL(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
const int *vartype=INTEGER(R_vartype);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%d ",vartype[i]);} Rprintf("\n");*/
network dag;
cycle cyclestore;
storage nodescore;
gsl_matrix *posterior;
/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
maxparents=asInteger(R_maxparents);
verbose=asInteger(R_verbose);
numvariates=asInteger(R_numvariates);/** number of points at which to evaluate the posterior distibution **/
whichnode=asInteger(R_whichnode); 
whichvariable=asInteger(R_whichvariable);
whichgaus=asInteger(R_whichgaus);
const int maxiters=asInteger(R_maxiters);
const double epsabs=asReal(R_epsabs);

/*Rprintf("got %d %d %d\n",whichnode,whichvariable, whichgaus);*/
SEXP listresults;
SEXP tmplistentry;
/** end of argument parsing **/

/*Rprintf("%f %f %f %f %f\n",priormean,priorsd,lowerlimit,upperlimit,relerr);*/
/** how to extract string entries of vectors of strings within a list */
/**Rprintf("->%s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),9))); the first entry in the list(), 
                                                                       then the 10 entry in that vector of strings 
                                                                       and this entry is a string**/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
PROTECT(listresults = allocVector(VECSXP, 2));
for(i=0;i<2;i++){
				PROTECT(tmplistentry=NEW_NUMERIC(numvariates));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm_mixed(R_obsdata,&obsdata, vartype);

/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag_mixed(&dag,&obsdata,maxparents);

/** all this does is to set internally set dag->defn[child][parent] etc **/
get_dag(&dag,R_dag);

/** check is the network definition is actually acyclic **/
init_hascycle(&cyclestore,&dag);

if(hascycle(&cyclestore,&dag)){error("network definition is not acyclic!\n");}

init_network_score(&nodescore,&dag);
/*calc_network_Score(&nodescore,&dag,&obsdata,priordatapernode, useK2,verbose, R_labels);*/
/** An R MATRIX it is single dimension and just needs to be unrolled */
posterior=gsl_matrix_alloc(numvariates,2);

gsl_set_error_handler_off();/*Rprintf("Note: turning off GSL Error handler\n");*/

for(j=0;j<2;j++){for(i=0;i<numvariates;i++){gsl_matrix_set(posterior,i,j,REAL(R_posterior)[i+j*numvariates]);}} 

/*for(i=0;i<numvariates;i++){Rprintf("%f %f\n",gsl_matrix_get(posterior,i,0),gsl_matrix_get(posterior,i,1));} */

calc_network_Marginals_laplace(&nodescore,&dag,&obsdata,verbose,&designmatrix,priormean,priorsd,priorgamshape,priorgamscale,
                               whichnode,whichvariable,whichgaus,posterior,maxiters,epsabs);
                                 
/** Roll back into an R MATRIX */
for(j=0;j<2;j++){for(i=0;i<numvariates;i++){REAL(VECTOR_ELT(listresults,j))[i]=gsl_matrix_get(posterior,i,j);}} 

gsl_matrix_free(posterior);

/*Rprintf("\n----------------------------------------------------------\n");
Rprintf("(LOG) NETWORK SCORE = %f\n",dag.networkScore);
Rprintf("----------------------------------------------------------\n");
*/


/*for(i=0;i<8;i++){x[i]= -x[i];}*/           
/*R_isort(x,8);
for(i=(len-1);i>=(len-5)+1-1;i--){Rprintf("%d\n",x[i]);}           
  */      
        /** set up memory storage for any network which has each node with <=maxparents **/

gsl_set_error_handler (NULL);/** restore the error handler*/

UNPROTECT(1);
/*PutRNGstate();*/

/*free_dag(&dag);*/
return(listresults);

}


