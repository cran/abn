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
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>

#define DEBUG_12

SEXP allnodesbn_additive(SEXP R_obsdata, SEXP R_priors_mean, SEXP R_priors_sd,SEXP R_priors_gamshape,SEXP R_priors_gamscale,
			 SEXP R_maxparents, SEXP R_verbose, SEXP R_vartype, SEXP R_maxiters, SEXP R_epsabs, SEXP R_errorverbose, SEXP R_whichnodes)
{
/** ****************/
/** declarations **/
int numNodes=LENGTH(R_obsdata);
unsigned int i,j,k,maxparents;
unsigned int verbose;
int errverbose;
datamatrix obsdata, designmatrix;
const double *priormean=REAL(R_priors_mean);/*Rprintf("priormean=%f %f\n",priormean[0],priormean[5]);*/
const double *priorsd=REAL(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double *priorgamshape=REAL(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double *priorgamscale=REAL(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
const int *vartype=INTEGER(R_vartype);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%u ",vartype[i]);} Rprintf("\n");*/
network dag;
storage nodescore;
/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
maxparents=asInteger(R_maxparents);
verbose=asInteger(R_verbose);
errverbose=asInteger(R_errorverbose);
SEXP listresults;
SEXP nodeidvector,scorevector,parentsvector, maxparentsvector;
const int maxiters=asInteger(R_maxiters);
const double epsabs=asReal(R_epsabs);
int curnode=0;
int curparentindex=0;
gsl_combination *c;
int total=0;
double totaldbl=0.0;
double curscore=0.0;
int *ptr_parents,*ptr_nodeid;
double *ptr_scores;
int localindex=0;/** this indexes every individual entry (not just row) for all parent combinations */
/** end of argument parsing **/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** need total number of entries **/
totaldbl=0.0;
for(i=0;i<=maxparents;i++){totaldbl+=gsl_sf_choose(numNodes-1,i);} /** total number of subsets per variable */
totaldbl= totaldbl*LENGTH(R_whichnodes);/** since we have total combinations per node */ 
/*Rprintf("total=%f %d %d %d\n",totaldbl,maxparents,numNodes,LENGTH(R_whichnodes));*/
total=rint(totaldbl);/** this is NEEDED - as gsl_sf_choose() returns double and just casting to int does not work - e.g. may be out by 1 rint() seems to work fine **/
/*Rprintf("total=%d %10.10f %d\n",total,totaldbl, (size_t)totaldbl);*/
PROTECT(listresults = allocVector(VECSXP, 4));/** number of elements in the outer most list**/
         PROTECT(nodeidvector=NEW_INTEGER(total));/** a single vector containing the node id from 1-n */
	 SET_VECTOR_ELT(listresults, 0, nodeidvector);
         PROTECT(scorevector=NEW_NUMERIC(total));/** a single vector containing the network score for each node and parent combination */
	 SET_VECTOR_ELT(listresults, 2, scorevector);/** assign the above vector to the first entry in the R list */
	 PROTECT(parentsvector = NEW_INTEGER(total*numNodes));/** one long vector since the matrix is filled by cols not rows anyway so have to transform (back in R) */
	 SET_VECTOR_ELT(listresults, 1, parentsvector);
	 PROTECT(maxparentsvector = NEW_INTEGER(1));/** just to hold maxparents*/
	 SET_VECTOR_ELT(listresults, 3, maxparentsvector);
	 UNPROTECT(4);
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */
gsl_set_error_handler_off();

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm_mixed(R_obsdata,&obsdata, vartype);
/** checked 21/05 - seems to work fine */
/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag_mixed(&dag,&obsdata,maxparents);
/** all this does is to set internally set dag->defn[child][parent] etc **/
init_network_score_mixed(&nodescore,&dag);

/** calc_network_Score_laplace(&nodescore,&dag,&obsdata,verbose,&designmatrix, priormean,priorsd,priorgamshape,priorgamscale,maxiters,epsabs,errverbose); **/

total=0; 
ptr_parents = INTEGER(parentsvector);/** point to the location of this matrix - actually one long vector **/
ptr_nodeid = INTEGER(nodeidvector);
ptr_scores = REAL(scorevector);
if(maxparents==(numNodes-1)){INTEGER(maxparentsvector)[0]=0;/** flag to say we are NOT using a parent restriction */
} else {INTEGER(maxparentsvector)[0]=1;} /** using a parent restriction **/

/** calc all possible parent combinations **/
for(i=0;i<LENGTH(R_whichnodes);i++){/** for each node passed get all network scores for all possible subsets subject to cardinality maxparents**/
  Rprintf("processing node...%d\n",INTEGER(R_whichnodes)[i]);
  curnode=INTEGER(R_whichnodes)[i]-1;/** -1 since R indexes start at unity and not zero **/ 
  /** get all possible subsets **/

  for(j=0; j<=maxparents; j++){
    c = gsl_combination_calloc (numNodes, j);/** setup for all combinations of (maxparents choose N) 
                                                 NOTE: use N and not N-1 since it means not having to remap the indices but we then need to discard
                                                 any combination where the parent=1 and the parent=curnode **/
    do
       {
       for(k=0;k<numNodes;k++){dag.defn[curnode][k]=0;}/** reset the node parents to independence **/
	 /** set the parents to the new combination **/
       for(k=0;k<gsl_combination_k(c);k++){curparentindex=gsl_combination_get (c, k);
	                                   dag.defn[curnode][curparentindex]=1;/** set parent */
                                          }
       if(dag.defn[curnode][curnode]!=0){continue;} /** this combination has the child as its own parent to skip to next iteration */
	 
       for(k=0;k<numNodes;k++){ptr_parents[localindex++]=dag.defn[curnode][k];}/** copy to store **/
                                                                  
       ptr_nodeid[total]=curnode+1;/** note the +1: back to R numbering */
       /*Rprintf("total STORE=%d\n",total);*/
      /** now if the current combination has more parents than maxparents then set score to a dummy value of zero (e.g. exp(0)=1) */
 
	switch(obsdata.vartype[curnode])  /** choose which type of node we have */
                 {
                      case 1:{ /** binary/categorical node */
                             curscore=calc_node_Score_laplace(&nodescore,&dag,&obsdata,curnode,errverbose, &designmatrix, priormean, priorsd,
                                                                  priorgamshape,priorgamscale,maxiters,epsabs); 
                             if(verbose){Rprintf("Binary node=%d score=%f\n", curnode,curscore);}                                                                                         
                             break;
                              }
                         
                      case 0:{ /** gaussian node */
                              /** note gaussiannodeid=0 - e.g. this uses the FIRST entry in each of the prior vectors for all gaussian nodes */
                              curscore=calc_Gaussiannode_Score_laplace(&nodescore,&dag,&obsdata,curnode,errverbose, &designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,0,maxiters,epsabs);
                              if(verbose){Rprintf("Gaussian node=%d score=%f\n",curnode,curscore);}
                              break;
                              }
                          default: {error("in default switch - should never get here!");}                                          
                   }                                     
                         
      ptr_scores[total]=curscore;/*if(curscore== -DBL_MAX){Rprintf("got dbl_max %d\n",total);}*/
       
      if(verbose){Rprintf("%d ",total);for(k=0;k<numNodes;k++){Rprintf("%d,",dag.defn[curnode][k]);} Rprintf("%10.10f\n",curscore);}
      total++;
      /** now calculate the network score **/
      
      /** we now have a new set of parents for node curnode so find score */
       }
       while (gsl_combination_next (c) == GSL_SUCCESS);
       gsl_combination_free (c);
    
       R_CheckUserInterrupt();/** allow an interupt from R console */
    } 
    
 } /** end of for each node passed */


gsl_set_error_handler (NULL);/** restore the error handler*/
UNPROTECT(1);
return(listresults);

}


