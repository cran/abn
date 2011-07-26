#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "utility_fns.h"
#include "network.h"

#include "laplace.h"
#include "laplace_marginals.h"

#define DEBUG_12

SEXP fitnetwork_additive(SEXP R_obsdata, SEXP R_dag,SEXP R_priors_mean, SEXP R_priors_sd,
			 SEXP R_maxparents,SEXP R_numVarLevels, SEXP R_labels, SEXP R_verbose)
{
/** ****************/
/** declarations **/
unsigned int i,maxparents;
unsigned int verbose;
datamatrix obsdata, designmatrix;
const double *priormean=REAL(R_priors_mean);/*Rprintf("priormean=%f\n",priormean[0]);*/
const double *priorsd=REAL(R_priors_sd);/*Rprintf("priorsd=%f\n",priorsd[0]);*/
network dag;
cycle cyclestore;
storage nodescore;
/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
maxparents=asInteger(R_maxparents);
verbose=asInteger(R_verbose);
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
PROTECT(listresults = allocVector(VECSXP, 1));
for(i=0;i<1;i++){
				PROTECT(tmplistentry=NEW_NUMERIC(1));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm(R_obsdata,&obsdata, R_numVarLevels);
/** checked 21/05 - seems to work fine */

/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag(&dag,&obsdata,maxparents);

/** all this does is to set internally set dag->defn[child][parent] etc **/
get_dag(&dag,R_dag);

/** check is the network definition is actually acyclic **/
init_hascycle(&cyclestore,&dag);

if(hascycle(&cyclestore,&dag)){error("network definition is not acyclic!\n");}

init_network_score(&nodescore,&dag);

/*calc_network_Score(&nodescore,&dag,&obsdata,priordatapernode, useK2,verbose, R_labels);*/

calc_network_Score_laplace(&nodescore,&dag,&obsdata,verbose,&designmatrix, R_labels,priormean,priorsd);

/*Rprintf("\n----------------------------------------------------------\n");
Rprintf("(LOG) NETWORK SCORE = %f\n",dag.networkScore);
Rprintf("----------------------------------------------------------\n");
*/
REAL(VECTOR_ELT(listresults,0))[0]=dag.networkScore;


/*for(i=0;i<8;i++){x[i]= -x[i];}*/           
/*R_isort(x,8);
for(i=(len-1);i>=(len-5)+1-1;i--){Rprintf("%d\n",x[i]);}           
  */      
        /** set up memory storage for any network which has each node with <=maxparents **/


UNPROTECT(1);
/*PutRNGstate();*/

/*free_dag(&dag);*/
return(listresults);

}


