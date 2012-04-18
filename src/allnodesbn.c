#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "utility_fns.h"
#include "network.h"
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>

#define DEBUG_12

SEXP allnodesbn(SEXP R_obsdata, SEXP R_useK2,SEXP R_maxparents,SEXP R_priorpernode, SEXP R_numVarLevels, SEXP R_verbose, SEXP R_whichnodes)
{
  
/** ****************/
/** declarations **/
unsigned int /*numObs,numNodes,*/i,j,k,maxparents;
unsigned int useK2, verbose;
double priordatapernode;
datamatrix obsdata;
network dag;
storage nodescore;
/** end of declarations*/
/** *******************/
/*GetRNGstate();*/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
int numNodes=LENGTH(R_obsdata);/*Rprintf("number of nodes=%d\n",numNodes);*/
/*numObs=LENGTH(VECTOR_ELT(R_obsdata,0));
obsdata.numVars=numNodes;*/
maxparents=asInteger(R_maxparents);
useK2=asInteger(R_useK2);
verbose=asInteger(R_verbose);
priordatapernode=asReal(R_priorpernode);
SEXP listresults;
SEXP nodeidvector,scorevector,parentsvector, maxparentsvector;
SEXP R_labels=0;/** note: this is no longer used but here as a dummy to keep other function calls happy**/
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

/** how to extract string entries of vectors of strings within a list */
/**Rprintf("->%s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),9))); the first entry in the list(), 
                                                                       then the 10 entry in that vector of strings 
                                                                       and this entry is a string**/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** need total number of entries **/
totaldbl=0.0;
for(i=0;i<=maxparents;i++){totaldbl+=gsl_sf_choose(numNodes-1,i);} /** total number of subsets per variable */
totaldbl= totaldbl*LENGTH(R_whichnodes);/** since we have total combinations per node */ 
total=rint(totaldbl);/** this is NEEDED - as gsl_sf_choose() returns double and just casting to int does not work - e.g. may be out by 1 rint() seems to work fine **/
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

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm(R_obsdata,&obsdata, R_numVarLevels);
/** checked 21/05 - seems to work fine */

/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag(&dag,&obsdata,maxparents);
init_network_score(&nodescore,&dag);
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

  for(j=0; j<=maxparents; j++){/**  */
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
      /** now if the current combination has more parents than maxparents then set score to a dummy value of zero (e.g. exp(0)=1) */
      curscore=calc_node_Score(&nodescore,&dag,&obsdata,curnode,priordatapernode, useK2, 0, R_labels);
              ptr_scores[total]=curscore;
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
  
	 
UNPROTECT(1);
return(listresults);

}


