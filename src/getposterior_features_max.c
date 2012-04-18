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

#define PRINTMEe
#define NotLOGa
#define UseLOG

SEXP getposterior_features_max(SEXP R_localscoreslist, SEXP R_numnodes, SEXP R_maxparents, SEXP R_offset)
{
/** ****************/
/** declarations **/
unsigned int numRows,i,j,index,k,l,found;
unsigned int numNodes,maxparents;
double *f_hat_0,*f_hat_1;
/** end of declarations*/
/** ********************************************/
/** parse function arguments - R data structs **/
numRows=LENGTH(VECTOR_ELT(R_localscoreslist,0));/** number of local node scores */ /*Rprintf("no. of score=%d\n",numRows);*/
numNodes=asInteger(R_numnodes);/** number of nodes in the network */
maxparents=asInteger(R_maxparents);
SEXP listresults,tmplistentry;
double  offset=asReal(R_offset);/** used to help avoid overflows **/
/** end of argument parsing **/

int *ptr_nodeid=INTEGER(VECTOR_ELT(R_localscoreslist,0));/** vector of node ids NOTE from 1 not zero **/
double *ptr_score=REAL(VECTOR_ELT(R_localscoreslist,2)); /** vector of node scores */

int **parents,*parentstmp;
double tmp=0.0;
#ifdef NotLOG
gsl_sf_result res;
int overflow=0;
#endif 
int numsubsets=0;
int **parents_short;/* **parents_short_ptr;*/
double **alpha,*alphatmp;
int **alphalookup,*alphalookuptmp;
int *V,lenV,toplevellen,*V2, *Vcopy;
int *parents_numparents,bestparents;
int tot=0;
int maxnode,indextmp;
int NumIter,TotalIters;
/*Rprintf("note: using logs - experimental...\n");*/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
PROTECT(listresults = allocVector(VECSXP, numNodes));
for(i=0;i<numNodes;i++){
				PROTECT(tmplistentry=NEW_INTEGER(numNodes+1));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
/** *******************************************************************************/

/** copy the R matrix - which is really a long vector - into a new 2-d array as it will produce easier to read code */
parents=(int **)R_alloc( (numRows),sizeof(int*));/** number of ROWS*/
	for(i=0;i<numRows;i++){parentstmp=(int *)R_alloc( numNodes,sizeof(int)); parents[i]=parentstmp;} 
/** now unroll the R matrix into the 2d array */
for(i=0;i<numRows;i++){for(j=0;j<numNodes;j++){parents[i][j]=INTEGER(VECTOR_ELT(R_localscoreslist,1))[i+numRows*j];}}
  

/** for each parent combination we need a prior value - q^\prime in the koivisto paper. Use the same prior as their where the prob of each parent combination
 *  for each node is proportional to 1/(choose(n-1,cardinality of G), e.g. G is the parent combination. Iterate over each entry in ptr_score[] and "update" this
  with prior*score*f(g) */

gsl_set_error_handler_off();/*Rprintf("Note: turning off GSL Error handler\n");*//** needed in case of underflow - use a manual catch */

for(i=0;i<numRows;i++){/** for every combination of parents at each node (subject to a cardinality constraint) */
  
   /*Rprintf("%d %d %d %d %d %d %5.5f %f %f\n",ptr_nodeid[i],parents[i][0],parents[i][1],parents[i][2],parents[i][3],parents[i][4],ptr_score[i],
	   qprime(parents[i],numNodes,qprimedenominator),feature(parents[i],ptr_nodeid[i]-1,child-1,parent-1));*/
   /** calculate the beta_i(G_i) values p557 in Koivisto - needed for each part of the later alpha calc**/
   tmp= qprime(parents[i],numNodes,offset);/** unlike in getposterior_features() we now have feature=1 since consider all features  */
   
#ifdef NotLOG
   overflow=gsl_sf_exp_mult_e(ptr_score[i],tmp, &res);/** raise ptr_score[i] to exp() then multiply by tmp; Overwrite back into score vector*/
   if(overflow){error("Error: %s\n",gsl_strerror(overflow));}
   ptr_score[i]= res.val;/** NOTE: using logs rather than transforming back to probs - since maximization **/
#endif   

#ifdef UseLOG
ptr_score[i]+= gsl_sf_log(tmp);
#endif

  /* Rprintf("Beta_i(G_i)=%20.20e\n",ptr_score[i]);*/
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

/** need an alpha store e.g. alpha[i][] is alpha_i and alpha[i][j] is the jth member of alpha_i **/
alpha=(double **)R_alloc( (numNodes),sizeof(double*));/** one for each node*/
	for(i=0;i<numNodes;i++){alphatmp=(double *)R_alloc( numsubsets,sizeof(double)); alpha[i]=alphatmp;} /** one for each subset **/
	  
f_hat_0=(double *)R_alloc( numsubsets,sizeof(double));/** enough space for each single alpha_i(.) */
f_hat_1=(double *)R_alloc( numsubsets,sizeof(double));/** enough space for each single alpha_i(.) */

for(i=0;i<numNodes;i++){
/** the ptr_score[] and ptr_nodeid[] need changed each time - results in f_hat_0[] 
    note - this returns the MAX value over each set of subsets and NOT the sum. Then need to find which subset corresponds to that sum  **/
mobius_transform_max(f_hat_0,f_hat_1,parents_short,ptr_score+numsubsets*i, numsubsets,ptr_nodeid+numsubsets*i, numNodes,alpha);/** pass a pointer to the start of the next node, f_hat_1 will have the results in **/
}
/** now for each entry in alpha_i[j] need to find which subset that corresponds to e.g. which parent pattern  **/
/** easiest to do this by creating a similar array for alpha[][] but where this now holds the row in parents of the relevant parent combination */
/** need an alpha store e.g. alpha[i][] is alpha_i and alpha[i][j] is the jth member of alpha_i **/
alphalookup=(int **)R_alloc( (numNodes),sizeof(int*));/** one for each node*/
	for(i=0;i<numNodes;i++){alphalookuptmp=(int *)R_alloc( numsubsets,sizeof(int)); alphalookup[i]=alphalookuptmp;} /** one for each subset **/

R_CheckUserInterrupt();/** allow an interupt from R console */

#ifdef PRINTME
/** this part just does some printing **/
for(k=0;k<numNodes;k++){/** for each node - child - in alpha **/
for(i=0;i<numsubsets;i++){/** for each subset of parents **/
	/** now find the correct parent set for each **/
      found=0;	
      for(l=0;l<numRows;l++){/** search each combination until find a hit **/
	if(ptr_score[l]==alpha[k][i]){Rprintf(" alpha[child][parentset]=%5.10e %5.10e\t",alpha[k][i],ptr_score[l]);found=1;
	  for(j=0;j<numNodes;j++){Rprintf("%d ",parents[l][j]);}Rprintf(" nodeid=%d\n",ptr_nodeid[l]); 
	break;}
      }
      if(found==0){error("no match found!!\n");}
	}
}
#endif

/** now need to find which parent set corresponds to the alpha values and store this index (in parents[][]) in alphalookup. Need this to be able to
 *  recover which parents are the optimal  ********************************************************************************************************/

for(k=0;k<numNodes;k++){/** for each node - child - in alpha **/
for(i=0;i<numsubsets;i++){/** for each subset of parents **/
	/** now find the correct parent set for each **/
      found=0;	
      for(l=0;l<numRows;l++){
	if(ptr_score[l]==alpha[k][i]){alphalookup[k][i]=l;found=1;break;} /** found a match so skip to next subset **/
      }
      if(found==0){error("no match found!!\n");}
	}
}

#ifdef PRINTME
for(k=0;k<numNodes;k++){/** for each node - child - in alpha **/
for(i=0;i<numsubsets;i++){/** for each subset of parents **/
  Rprintf(" alphalookup[child][parentset]=%d\n",alphalookup[k][i]);}}
#endif

/** *************************************************************************************************************************************/
/** **** We now have alpha_i(S) for every i and S where this contains the maximum value and a lookup of the relevant parent combination */
/** **** in alpha[i][s] and alphalookup[i][s] respectively ******************************************************************************/
/** alpha is the form alpha[i][j] is the ith node (indexing from 0) and then the jth parent combination (indexing from 0) for that node */
/** We now find g(V) but where we want to find the MAX of the sum at each stage in the recursive g() calc. this is a pain but seems to **/
/** work, however it is probably inefficient. The hardest part is finding out which parents are those for each node and this is current */
/** implemented by re-running the recursion for successively smaller sets and each times gives us the optimal combination for one node  */
/** *************************************************************************************************************************************/
/** FIRST - always have one pass**/

TotalIters=numNodes-1; /** the last node is ALWAYS independent and so don't need an iteration for this **/ 
NumIter=1;
lenV=numNodes;toplevellen=numNodes;
V=(int *)R_alloc(lenV,sizeof(int));
for(i=0;i<lenV;i++){V[i]=i;} /** set up an arry of nodeids **/

/** first run on the full node set - this will return the node at the outermost level in maxnode  e.g. g(0,1,2)=alpha_0(1,2)g(1,2) so we
    now know that 0 is the first node and its parents are whatever alphalookup[0][parentset=12]. Note that parentset here will always be
    the biggest possible set excluding the child since we are working outwards in*/
g_max(V, lenV,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,0,toplevellen,&maxnode,ptr_score);
/*Rprintf("g_max=%5.10e\n",a);
Rprintf("first got maxnode=%d len=%d\n",maxnode,lenV);*/
indextmp=0;
V2=(int *)R_alloc(lenV-1,sizeof(int));for(i=0;i<lenV;i++){if(V[i]!=maxnode){V2[indextmp++]=V[i];}} 
bestparents=get_alpha_parents(maxnode, V, V2, lenV, alphalookup, parents, parents_numparents, ptr_nodeid, numRows, numNodes, maxparents, ptr_score);

/* Rprintf("nodeid=%d\t",ptr_nodeid[bestparents]);for(i=0;i<numNodes;i++){Rprintf("%d ",parents[bestparents][i]);}Rprintf("\n"); */

INTEGER(VECTOR_ELT(listresults,NumIter-1))[0]=ptr_nodeid[bestparents];/** store the child NODEID **/
for(i=0;i<numNodes;i++){INTEGER(VECTOR_ELT(listresults,NumIter-1))[i+1]=parents[bestparents][i];} /** store the parents **/

Vcopy=(int *)R_alloc(lenV,sizeof(int));for(i=0;i<lenV;i++){Vcopy[i]=V[i];}
/* idea is to re-run but dropped out the a variable each time since we know which is best */

while(NumIter<=TotalIters){
  
  lenV--;toplevellen--;
  V=(int *)R_alloc(lenV,sizeof(int));indextmp=0;
  for(i=0;i<lenV+1;i++){if(Vcopy[i]!=maxnode){V[indextmp++]=Vcopy[i];}}
  g_max(V, lenV,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,0,toplevellen,&maxnode, ptr_score);
  /*Rprintf("second got maxnode=%d len=%d\n",maxnode,lenV);*/
  indextmp=0;
  V2=(int *)R_alloc(lenV-1,sizeof(int));for(i=0;i<lenV;i++){if(V[i]!=maxnode){V2[indextmp++]=V[i];}} 
  bestparents=get_alpha_parents(maxnode, V, V2, lenV, alphalookup, parents, parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,ptr_score);
 /* Rprintf("nodeid=%d\t",ptr_nodeid[bestparents]);for(i=0;i<numNodes;i++){Rprintf("%d ",parents[bestparents][i]);}Rprintf("\n"); */
  INTEGER(VECTOR_ELT(listresults,NumIter))[0]=ptr_nodeid[bestparents];
  for(i=0;i<numNodes;i++){INTEGER(VECTOR_ELT(listresults,NumIter))[i+1]=parents[bestparents][i];} /** store the parents **/
  Vcopy=(int *)R_alloc(lenV,sizeof(int));for(i=0;i<lenV;i++){Vcopy[i]=V[i];}
  
NumIter++;  
R_CheckUserInterrupt();/** allow an interupt from R console */
}

#ifdef JUNK
/** SECOND **/
lenV--;toplevellen--;
V=(int *)R_alloc(lenV,sizeof(int));indextmp=0;
for(i=0;i<lenV+1;i++){if(Vcopy[i]!=maxnode){V[indextmp++]=Vcopy[i];}}
/*for(i=0;i<lenV;i++){Rprintf("%d ",V[i]);}Rprintf("\n");*/

g_max(V, lenV,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,0,toplevellen,&maxnode);
Rprintf("second got maxnode=%d len=%d\n",maxnode,lenV);

indextmp=0;
V2=(int *)R_alloc(lenV-1,sizeof(int));for(i=0;i<lenV;i++){if(V[i]!=maxnode){V2[indextmp++]=V[i];}} 
bestparents=get_alpha_parents(maxnode, V, V2, lenV, alphalookup, parents, parents_numparents, ptr_nodeid, numRows, numNodes, maxparents);

Rprintf("nodeid=%d\t",ptr_nodeid[bestparents]);for(i=0;i<numNodes;i++){Rprintf("%d ",parents[bestparents][i]);}Rprintf("\n");

Vcopy=(int *)R_alloc(lenV,sizeof(int));for(i=0;i<lenV;i++){Vcopy[i]=V[i];}


lenV--;toplevellen--;
V=(int *)R_alloc(lenV,sizeof(int));indextmp=0;
for(i=0;i<lenV+1;i++){if(Vcopy[i]!=maxnode){V[indextmp++]=Vcopy[i];}}
/*for(i=0;i<lenV;i++){Rprintf("%d ",V[i]);}Rprintf("\n");*/

g_max(V, lenV,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,0,toplevellen,&maxnode);
Rprintf("third got maxnode=%d len=%d\n",maxnode,lenV);

indextmp=0;
V2=(int *)R_alloc(lenV-1,sizeof(int));for(i=0;i<lenV;i++){if(V[i]!=maxnode){V2[indextmp++]=V[i];}} 
bestparents=get_alpha_parents(maxnode, V, V2, lenV, alphalookup, parents, parents_numparents, ptr_nodeid, numRows, numNodes, maxparents);

Rprintf("nodeid=%d\t",ptr_nodeid[bestparents]);for(i=0;i<numNodes;i++){Rprintf("%d ",parents[bestparents][i]);}Rprintf("\n");

Vcopy=(int *)R_alloc(lenV,sizeof(int));for(i=0;i<lenV;i++){Vcopy[i]=V[i];}

#endif
/*Rprintf("leV=%d\n",lenV);*/
/*
V=(int *)R_alloc(lenV,sizeof(int));indextmp=0;
for(i=0;i<lenV+1;i++){if(Vcopy[i]!=maxnode){V[indextmp++]=i;}}
for(i=0;i<lenV;i++){Rprintf("%d ",V[i]);}Rprintf("\n");

g_max(V, lenV,alpha,parents,parents_numparents, ptr_nodeid, numRows, numNodes, maxparents,0,toplevellen,&maxnode);
Rprintf("third got maxnode=%d\n",maxnode);

indextmp=0;
V2=(int *)R_alloc(lenV-1,sizeof(int));for(i=0;i<lenV;i++){if(V[i]!=maxnode){V2[indextmp++]=V[i];}} 
bestparents=get_alpha_parents(maxnode, V, V2, lenV, alphalookup, parents, parents_numparents, ptr_nodeid, numRows, numNodes, maxparents);

Rprintf("nodeid=%d\t",ptr_nodeid[bestparents]);for(i=0;i<numNodes;i++){Rprintf("%d ",parents[bestparents][i]);}Rprintf("\n");

Vcopy=(int *)R_alloc(lenV,sizeof(int));for(i=0;i<lenV;i++){Vcopy[i]=V[i];}
*/
  /*Rprintf("feature posterior prob=%20.20e\n",feature_postprob);*/

gsl_set_error_handler (NULL);/** restore the error handler*/   

UNPROTECT(1);

return(listresults);

}


