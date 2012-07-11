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
#include "mobius.h"

#define DEBUG_12

SEXP getposterior_features(SEXP R_localscoreslist, SEXP R_numnodes,SEXP R_childnode,SEXP R_parentnode, SEXP R_maxparents, SEXP R_featuredefn, SEXP R_offset)
{
/** ****************/
/** declarations **/
unsigned int numRows,i,j,index;
unsigned int child, parent, numNodes,maxparents;
double *f_hat_0,*f_hat_1;
/** end of declarations*/
/** ********************************************/
/** parse function arguments - R data structs **/
numRows=LENGTH(VECTOR_ELT(R_localscoreslist,0));/** number of local node scores */ /*Rprintf("no. of score=%d\n",numRows);*/
child=asInteger(R_childnode);
parent=asInteger(R_parentnode);
numNodes=asInteger(R_numnodes);/** number of nodes in the network */
maxparents=asInteger(R_maxparents);
SEXP listresults,tmplistentry;
double  offset=asReal(R_offset);/** used to help avoid overflows **/
/** end of argument parsing **/

int *ptr_nodeid=INTEGER(VECTOR_ELT(R_localscoreslist,0));/** vector of node ids NOTE from 1 not zero **/
double *ptr_score=REAL(VECTOR_ELT(R_localscoreslist,2)); /** vector of node scores */
int *customfeature=0;
if(!isNull(R_featuredefn)){customfeature=INTEGER(R_featuredefn);}/** this will hold a vector with a 1 or 0 for each and every possible parent combination **/

int **parents,*parentstmp;
double tmp=0.0;
gsl_sf_result res;
int overflow=0;
int numsubsets=0;
int **parents_short;/* **parents_short_ptr;*/
double **alpha,*alphatmp;
int *V;
int *parents_numparents;
int tot=0;
double feature_postprob;

/*gsl_sf_result_e10 bignum;
double x= -10.0;
double y= 3.0;
gsl_sf_exp_mult_e10_e (x,y, &bignum);

Rprintf("==%5.10e %5.10e\n",pow(10.0,0.0001361),pow(9.5,3));
Rprintf("%5.10e %d %5.10e\n",bignum.val,bignum.e10, bignum.val*pow(10.0,(double)bignum.e10));
*/
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
/** *******************************************************************************/
/*Rprintf("%d %d %f\n",INTEGER(VECTOR_ELT(R_localscoreslist,0))[0],INTEGER(VECTOR_ELT(R_localscoreslist,0))[16],REAL(VECTOR_ELT(R_localscoreslist,2))[0]);*/

/** copy the R matrix - which is really a long vector - into a new 2-d array as it will produce easier to read code */
parents=(int **)R_alloc( (numRows),sizeof(int*));/** number of ROWS*/
	for(i=0;i<numRows;i++){parentstmp=(int *)R_alloc( numNodes,sizeof(int)); parents[i]=parentstmp;} 
/** now unroll the R matrix into the 2d array */
for(i=0;i<numRows;i++){for(j=0;j<numNodes;j++){parents[i][j]=INTEGER(VECTOR_ELT(R_localscoreslist,1))[i+numRows*j];}}
  /*for(i=0;i<numRows;i++){for(j=0;j<numNodes;j++){Rprintf("%d ",parents[i][j]);}Rprintf("nodeid=%d score=%f\n",ptr_nodeid[i],ptr_score[i]);}*/ 

/** for each parent combination we need a prior value - q^\prime in the koivisto paper. Use the same prior as their where the prob of each parent combination
 *  for each node is proportional to 1/(choose(n-1,cardinality of G), e.g. G is the parent combination. Iterate over each entry in ptr_score[] and "update" this
  with prior*score*f(g) */

/*qprimedenominator=0;
for(i=0;i<=maxparents; i++){qprimedenominator+= gsl_sf_choose(numNodes-1,i);}
*/
gsl_set_error_handler_off();/*Rprintf("Note: turning off GSL Error handler\n");*//** needed in case of underflow - use a manual catch */

/*Rprintf("\n ---- NOTE - still using dummy values for BETA ----\n");*/

for(i=0;i<numRows;i++){/** for every combination of parents at each node (subject to a cardinality constraint) */
  
   /*Rprintf("%d %d %d %d %d %d %5.5f %f %f\n",ptr_nodeid[i],parents[i][0],parents[i][1],parents[i][2],parents[i][3],parents[i][4],ptr_score[i],
	   qprime(parents[i],numNodes,qprimedenominator),feature(parents[i],ptr_nodeid[i]-1,child-1,parent-1));*/
   /** calculate the beta_i(G_i) values p557 in Koivisto - needed for each part of the later alpha calc**/
   tmp= qprime(parents[i],numNodes,offset)*feature(parents[i],ptr_nodeid[i],child,parent,i,customfeature,numRows);/**  */
   if(ptr_score[i]== -DBL_MAX){ptr_score[i]=0.0;Rprintf("Warning: got a -DBL_MAX\n");
   } else {
     overflow=gsl_sf_exp_mult_e(ptr_score[i],tmp, &res);
     /*overflow=gsl_sf_exp_mult_e10_e(ptr_score[i],tmp, &bignum);*//** raise ptr_score[i] to exp() then multiply by tmp; Overwrite back into score vector*/
     if(overflow){error("Error: %s\n",gsl_strerror(overflow));
                /*overflow=gsl_sf_exp_mult_e10_e(ptr_score[i],tmp, &bignum);*/    
  }	
   ptr_score[i]=res.val;
   /*ptr_score[i]=bignum.val* pow(10.0,(double)bignum.e10);*//** JUST STORE THE 10^ EXPONENT - as otherwise this will exceed floating point range! - a bit desperate ! */ 
   }
   
   /*Rprintf("Beta_i(G_i)=%20.20e\n",ptr_score[i]);*/ 
   }

/** we now have a vector ptr_score[] with all the values needed for the alpha_i(.) calc **/
/** now for fast mobius transform for EACH alpha_i(.) for i in 0 through n-1 **/ 
/** to make the next part easier to follow copy all the parent combinations into a new array with the column for the child node dropped in each case */
/** IMPORTANT NOTE: as the combinations for EACH node are identical once the child node has been removed only pass the combinations for a single NODE */
numsubsets=numRows/numNodes;
parents_short=(int **)R_alloc( (numsubsets),sizeof(int*));/** number of ROWS*/
	for(i=0;i<numsubsets;i++){parentstmp=(int *)R_alloc( numNodes-1,sizeof(int)); parents_short[i]=parentstmp;}
/*parents_short_ptr=parents_short;*//** store the address - needed later **/
for(i=0;i<numsubsets;i++){index=0;
                       for(j=0;j<numNodes;j++){if(j!=ptr_nodeid[i]){/** if NOT the child column */
			                                            parents_short[i][index++]=parents[i][j];}
	}
}

/** also need later the number of parents in each parent combination and do this now to avoid continually re-calculating this **/
parents_numparents=(int *)R_alloc( (numRows),sizeof(int));/** number of ROWS*/
for(i=0;i<numRows;i++){tot=0;for(j=0;j<numNodes;j++){if(parents[i][j]==1){tot++;}} parents_numparents[i]=tot;}

/*for(i=0;i<numsubsets;i++){for(j=0;j<numNodes-1;j++){Rprintf("%d ",parents_short[i][j]);}Rprintf(" numparent=%d\n",parents_numparents[i]);}  */

/** need an alpha store e.g. alpha[i][] is alpha_i and alpha[i][j] is the jth member of alpha_i **/
alpha=(double **)R_alloc( (numNodes),sizeof(double*));/** one for each node*/
	for(i=0;i<numNodes;i++){alphatmp=(double *)R_alloc( numsubsets,sizeof(double)); alpha[i]=alphatmp;} /** one for each subset **/
	  
f_hat_0=(double *)R_alloc( numsubsets,sizeof(double));/** enough space for each single alpha_i(.) */
f_hat_1=(double *)R_alloc( numsubsets,sizeof(double));/** enough space for each single alpha_i(.) */

for(i=0;i<numNodes;i++){
/** the ptr_score[] and ptr_nodeid[] need changed each time - results in f_hat_0[] **/
mobius_transform(f_hat_0,f_hat_1,parents_short,ptr_score+numsubsets*i, numsubsets,ptr_nodeid+numsubsets*i, numNodes,alpha);
/** pass a pointer to the start of the next node, f_hat_1 will have the results in **/
}
/*parents_short=parents_short_ptr;*//** reset the address back to start of array */

/*
for(k=0;k<numNodes;k++){
for(i=0;i<numsubsets;i++){Rprintf("node:%d| ",k);
                          for(j=0;j<numNodes-1;j++){
                          Rprintf("%d ",parents_short[i][j]);}
			  Rprintf(" alpha[node][subset]=%5.10e\n",alpha[k][i]);}
}
*/
/** alpha is the form alpha[i][j] is the ith node (indexing from 0) and then the jth parent combination (indexing from 0) for that node **/

V=(int *)R_alloc(numNodes,sizeof(int));
for(i=0;i<numNodes;i++){V[i]=i;} /** set up an arry of nodeids **/

/*alpha[0][0]=1.1; alpha[0][1]=0.3; alpha[0][2]=0.5;  alpha[0][3]=0.1;  
alpha[1][0]=1.2; alpha[1][1]=0.6; alpha[1][2]=0.7;  alpha[1][3]=0.2; 
alpha[2][0]=1.3; alpha[2][1]=0.8; alpha[2][2]=0.9;  alpha[2][3]=0.3; 
*/

feature_postprob=g(V, numNodes,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,ptr_score);
/*Rprintf("feature posterior prob=%20.20e\n",feature_postprob);*/

REAL(VECTOR_ELT(listresults,0))[0]=feature_postprob;

gsl_set_error_handler (NULL);/** restore the error handler*/   

UNPROTECT(1);

return(listresults);

}


