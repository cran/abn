#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include "mobius.h"

#define NotLOGa
#define UseLOG

double qprime(int *parents, int numNodes, double offset)
{
  /*Rprintf("=>%d %d %d %d %d\n",parents[0],parents[1],parents[2],parents[3],parents[4]);*/
 int numparents=0;
 int i;
 for(i=0;i<numNodes;i++){numparents+=parents[i];} /** get number of parents in this combination */
 /*Rprintf("num pars=%d\n",numparents);*/ 
 /** number of ways to choose a combination with numparents out of a total of (length-1 parents then divide by constant to get a prob.mass.fn**/
 
 return(offset*1.0/gsl_sf_choose(numNodes-1,numparents));
  
}

/** **********************************************************************************/
double feature(int *parents, int nodeid, int child, int parent, int i, int *customfeature, int numRows)
{

  if(child== -1){/** this is a flag for use with calculating the demoninator e.g. f_i(G_i)=1 for all i and all G_i **/
                 return(1.0);}

  if(customfeature!=0){return(1.0*customfeature[i]);}               
                 
  if(nodeid!=child){/** not the child and so include by default **/
                       return(1.0);
  } else { /** have correct child so find out if arc parent->child is present in this parent combination */
         if(parents[parent]==1){return(1.0);/** arc is present to include */
	} else {return(0.0);/** arc is not present so exclude */
	} 
  }
  
}
/** **********************************************************************************/
void mobius_transform(double *f_hat_0,double *f_hat_1, int **parentsloc, double *score, int numsubsets, int *nodeid, int numnodes, double **alpha){

int i,j;
int n=numnodes-1;
/** initial step - f_hat_0=ptr_score for all i **/
for(i=0;i<numsubsets;i++){f_hat_0[i]=score[i];}

/*Rprintf("passed: ");for(k=0;k<numsubsets;k++){Rprintf("%20.20e ",score[k]);} Rprintf("\n"); */ 

/** consider subsets of 1 through n */
for(j=1;j<=n;j++){/** for each iteration in the mobius transform */
  /*Rprintf("iteration %d\n",j);*/
  for(i=0;i<numsubsets;i++){/** for each subset */
    /*index=0;
    for(k=0;k<numnodes;k++){if(k!=nodeid){subset[index++]=parentsloc[i][k];}} *//** copy each set of parents into new array dropping the column for the child */
    /*Rprintf("passed: ");for(k=0;k<numnodes;k++){Rprintf("%d ",parentsloc[i][k]);} Rprintf("\n");*/  
     /* if(index!=numnodes-1){error("ERROR\n");} */ /** a simple check */
    if(parentsloc[i][j-1]==0){/** if "j" is not a member of the current subset**/
             f_hat_1[i]= f_hat_0[i];/** f_{n}=f_{n-1} */
    } else { /** if "j" is a member of the current subset then have two parts f_{n-1}(X\{j}) and f_{n-1}(X) **/
            f_hat_1[i]= f_hat_0[i] + f_hat_0[indexjcomplement(i,j-1,parentsloc, n,numsubsets)];  
    }
  }
  /** now have f_hat_1[] for all subsets so copy into f_hat_0[] and repeat **/ 
 for(i=0;i<numsubsets;i++){f_hat_0[i]=f_hat_1[i];}
  }

 /** copy into alpha[] **/
 for(i=0;i<numsubsets;i++){alpha[*nodeid][i]=f_hat_0[i];}
  
}

/** **********************************************************************************/
/** rather than find sum over subsets find subset with the max value **/
void mobius_transform_max(double *f_hat_0,double *f_hat_1, int **parentsloc, double *score, int numsubsets, int *nodeid, int numnodes, double **alpha){

int i,j;
int n=numnodes-1;
double tmp;
/** initial step - f_hat_0=ptr_score for all i **/
for(i=0;i<numsubsets;i++){f_hat_0[i]=score[i];}

/*Rprintf("passed: ");for(k=0;k<numsubsets;k++){Rprintf("%20.20e ",score[k]);} Rprintf("\n"); */ 

/** consider subsets of 1 through n */
for(j=1;j<=n;j++){/** for each iteration in the mobius transform */
  /*Rprintf("iteration %d\n",j);*/
  for(i=0;i<numsubsets;i++){/** for each subset */
    /*index=0;
    for(k=0;k<numnodes;k++){if(k!=nodeid){subset[index++]=parentsloc[i][k];}} *//** copy each set of parents into new array dropping the column for the child */
    /*Rprintf("passed: ");for(k=0;k<numnodes;k++){Rprintf("%d ",parentsloc[i][k]);} Rprintf("\n");*/  
     /* if(index!=numnodes-1){error("ERROR\n");} */ /** a simple check */
    if(parentsloc[i][j-1]==0){/** if "j" is not a member of the current subset**/
             f_hat_1[i]= f_hat_0[i];/** f_{n}=f_{n-1} */
    } else { /** if "j" is a member of the current subset then have two parts f_{n-1}(X\{j}) and f_{n-1}(X) **/
            /** the MAX PART -  do this completely analogous to the sum transform just replace "+" with "max" */
	    tmp=f_hat_0[indexjcomplement(i,j-1,parentsloc, n,numsubsets)];
	    if(f_hat_0[i]>tmp){f_hat_1[i]=f_hat_0[i];} else {f_hat_1[i]=tmp;} /** simply choose the max value of either f^hat_{j-1}(X\{j}) or f^hat_{j-1}(X) **/
    }
  }
  /** now have f_hat_1[] for all subsets so copy into f_hat_0[] and repeat **/ 
 for(i=0;i<numsubsets;i++){f_hat_0[i]=f_hat_1[i];}
  }

 /** copy into alpha[] **/
 for(i=0;i<numsubsets;i++){alpha[*nodeid][i]=f_hat_0[i];}
  
}
/** ************************************************************************************/
int indexjcomplement(int currow, int curparent, int **parents, int numnodes, int numsubsets)
{
  
/** find the index of the row in parents[][] which has the combination equal to parents[currow][] but with the "curparent"th entry 0 */
int i,j,match;

/*Rprintf("got: ");for(i=0;i<numnodes;i++){Rprintf("%d ",parents[currow][i]);}Rprintf("\n"); */
 
for(i=0;i<numsubsets;i++){/** for each possible parent combination **/
  if(parents[i][curparent]==0){/** parent of interest is not in this parent set so check rest of this combination **/  
  /*for(j=0;j<numnodes;j++){*//** so are the rest of the parents the same as the parent set of interest*/
  /*  Rprintf("%d ",parents[i][j]);
  } Rprintf("\n");*/
  match=0;
  for(j=0;j<numnodes;j++){if(j!=curparent){
                                           if(parents[i][j]==parents[currow][j]){match=0;} else {match=1;break;}
  }}
  /*Rprintf("match=%d\n",match);*/
  if(match==0){/*Rprintf("i=%d\n",i);*/return(i);}  
  }
}

/** if get to here we have a problem! **/
error("subset not found!");
return(-1);
}
/** ******************************************************************************************/
/** Define the recursion function g(S) in Koivisto and Sood - V is the array of all node ids */
/** note this is a bit tricky and has been checked manually against V={0,1,2} etc            */
/** ******************************************************************************************/
double g(int *V, int len, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, double *ptr_score){
int i,j;
int *V2;
int index=0;
double tot=0.0;
double alpha_val=0.0;
int found;

 if(len==0){return(1);} /** this is the boundary case - e.g. g(empty_set)=1 */
 
 for(i=0;i<len;i++){/** for each member of set V **/
                    /*V2=(int *)R_alloc(len-1,sizeof(int));*//** allocate new smaller array. n.b. this seems wasteful in terms of memory but R will reclaim this since using R_alloc() */
		    V2=(int *)malloc((len-1)*sizeof(int));
		    index=0;
		    for(j=0;j<len;j++){if(j!=i){V2[index++]=V[j];}} /** copy entries into the new array **/
		    if(len-1<=maxparents){/** e.g. no parent limit lenV-1 is the length of V2**/
		    alpha_val = get_alpha(i,V,V2,len,alpha,parents, parents_numparents, nodesid, numRows, numNodes, maxparents, &found);
		    } else {/** we are asking for alpha_i(S) where |S| is > maxparents so need to do something different since this doesn't use mobius */
		    alpha_val = get_alpha_no_mobius(V[i],V2,parents,parents_numparents,nodesid,ptr_score,numRows,numNodes,maxparents,len-1);}
		       tot+=get_q(i,V,V2,len)*
		          alpha_val*
		          g(V2,len-1,alpha,parents, parents_numparents, nodesid, numRows,numNodes, maxparents, ptr_score);/** q is prior over orders, alpha is result of the mobius transform **/
	            free(V2);	    
        }
   
return(tot);  
}

/** ******************************************************************************************/
/** as above but adapted to maximise rather than sum - seems tricky but seems to work with literally just replacing */
/** the sum with a max() so that anywhere a sum would have happened the max individual element is returned       */
/** ******************************************************************************************/
double g_max(int *V, int len, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents,int previousg, 
	     int toplevellen, int *maxnode, double *ptr_score){
int i,j;
int *V2;
int index=0;
previousg++;
gsl_vector *totstore;
double maxtot=0.0;
int maxindex;
int found;
double alpha_val=0.0;

 if(len==0){/*Rprintf("empty\n");*/return(1);} /** this is the boundary case - e.g. g(empty_set)=1 */
 
 totstore = gsl_vector_alloc (len);
 
 for(i=0;i<len;i++){/** for each member of set V **/
                    /*V2=(int *)R_alloc(len-1,sizeof(int));*//** allocate new smaller array. n.b. this seems wasteful in terms of memory but R will reclaim this since using R_alloc() */
		    V2=(int *)malloc((len-1)*sizeof(int));
		    index=0;
		    for(j=0;j<len;j++){if(j!=i){V2[index++]=V[j];}} /** copy entries into the new array **/
		    if(len-1<=maxparents){/** e.g. no parent limit lenV-1 is the length of V2**/
		     alpha_val = get_alpha(i,V,V2,len,alpha,parents, parents_numparents, nodesid, numRows, numNodes, maxparents, &found);  
		     } else {/** we are asking for alpha_i(S) where |S| is > maxparents so need to do something different since this doesn't use mobius */
		      alpha_val = get_alpha_no_mobius_max(V[i],V2,parents,parents_numparents,nodesid,ptr_score,numRows,numNodes,maxparents,len-1);} 
#ifdef NotLOG
		     gsl_vector_set(totstore,i,
		         alpha_val*
		          g_max(V2,len-1,alpha,parents, parents_numparents, nodesid, numRows,numNodes, maxparents,previousg,toplevellen, maxnode,ptr_score)
				   );
#endif
		     
#ifdef UseLOG

 gsl_vector_set(totstore,i,
		         alpha_val+
		          g_max(V2,len-1,alpha,parents, parents_numparents, nodesid, numRows,numNodes, maxparents,previousg,toplevellen, maxnode,ptr_score)
				   );		     
#endif		     
		     
                   free(V2);
                  }
 
 maxtot=gsl_vector_max(totstore); /** whenever there is a sum over values then return the maximum of the individual values **/
 maxindex=gsl_vector_max_index(totstore);
 
 /*Rprintf("contents of totstore[]\n");for(i=0;i<len;i++){Rprintf("%5.5e ",gsl_vector_get(totstore,i));}Rprintf("\n");*/
 /** this is crucial - not just printing.... */
 if(len==toplevellen){/*Rprintf("\n\n ---- toplevel ----\n");
                      for(i=0;i<len;i++){Rprintf("%5.5e ",gsl_vector_get(totstore,i));}Rprintf("\n");
                      Rprintf("node is %d\n",V[maxindex]); */
		      *maxnode=V[maxindex];}
 
 gsl_vector_free(totstore);

R_CheckUserInterrupt();/** allow an interupt from R console */

return(maxtot);  
}

/** ***********************************************************************************************************************/
/** ***********************************************************************************************************************/
/** the prior over orders **/
double get_q(int dropped, int *V, int *V2, int lenV){
  
/** just assume uniform for the moment which means that for every choice of order (e.g. subset) this is simply 1/n! **/
/** obvious when writing out all the possible orders for three nodes etc. It is not 1/(n-1)! since the node may be first
 *  in the order (e.g. only the empty set before it. NOTE: we do NOT return 1/n! since this may be very small and then
 *  values of alpha (needed in g()) will already be potential underflows so just return unity as all the probabilities
 *  require to be standardized anyway so its not necessary for use of a denominator for the priors                     */ 
  return(1.0);
  
}
/** ***********************************************************************************************************************/
/** ***********************************************************************************************************************/
/** the key part of g()   *************************************************************************************************/

double get_alpha(int dropped, int *V, int *V2, int lenV, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, int *found)
{
  /** alpha is the form alpha[i][j] is a double value for the ith node (indexing from 0) and then the jth parent combination
   * (indexing from 0) for that node **/
  /** the key part of this function is simply to find the correct entry in alpha which has been reqested by g() **/
  
  /** NOTE: we know the node as this is V[dropped] so can go right to the correct row in alpha **/ 
 int i,j;
 int check=0;
 *found=1;
 /*Rprintf("lenV=%d\n",lenV);*/
       /*if(lenV==0){Rprintf("empty\n");}*/
       /*for(i=0;i<lenV;i++){Rprintf("%d ",V[i]);}Rprintf("\n");*/ /*best line */
       /*for(i=0;i<lenV-1;i++){Rprintf("%d ",V2[i]);}Rprintf("\n");*/
    
 /*if(lenV-1>maxparents){*//** have asked for a parent combinaton with more parents than allowed e.g. alpha(0,1,2) when parent limit is 2 etc **/
  /* return(get_alpha_maxparent_exceeded(dropped,V, V2, lenV, alpha, parents, parents_numparents, nodesid, numRows, numNodes, maxparents));
                    }*/
 
 for(i=0;i<numRows;i++){/** for each possible subset (e.g. alpha) - just need to locate the correct entry in alpha */ 
   /*Rprintf("i=%d\n",i);*/
   /** case 1. empty set **/
   if(   parents_numparents[i]==(lenV-1) /** got right length of subset, -1 as this is analogous to V2 */
      && parents_numparents[i]==0 /** and the empty set **/
      && nodesid[i]==V[dropped] /** got correct node */
   ){ /** we are done **/
      /*Rprintf("alpha[ V[dropped] ][i%(numRows/numNodes)]=%f\n",alpha[ V[dropped] ][i%(numRows/numNodes)]);*/
      return(alpha[ V[dropped] ][i%(numRows/numNodes)]);}
      
   if(   parents_numparents[i]==(lenV-1) /** -1 as this is analogous to V2 */
      && parents_numparents[i]!=0 
      && nodesid[i]==V[dropped]
    ){ check=0;
       for(j=0;j<lenV-1;j++){/** for each member of V2 **/
	 /*Rprintf("j=%d\n",j);*/
	   
	   if(parents[i][ V2[j] ]==1){check++;}
       }
       /*Rprintf("check=%d\n",check);*/
       if(check==parents_numparents[i]){/** we have both the correct number of parents and the correct combination of parents */
	 /*Rprintf("||");Rprintf("node=%d ",nodesid[i]);for(j=0;j<numNodes;j++){Rprintf("%d ",parents[i][j]);}Rprintf("||\n");
	 Rprintf("alpha[ V[dropped] ][i%(numRows/numNodes)]=%f\n",alpha[ V[dropped] ][i%(numRows/numNodes)]);*/
	 return(alpha[ V[dropped] ][i%(numRows/numNodes)]);}
   }
 }
 
 /** if we get to here is it because the number of parents is less than numNodes-1 e.g. a restricted parent set **/
 /** hence some of the alpha_i(S) will not exist **/
 *found=0;/** flag to say did not find alpha!! **/  /** "should never get here! **/

 return(-1.0);  

  
}
/** **************************************************************************************************************/
/** like get_alpha but returns an index to the parents of alpha_i(S) where these are the max subset of S         */
/** **************************************************************************************************************/
int get_alpha_parents(int dropped, int *V, int *V2, int lenV, int **alphalookup, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, double *ptr_score)
{
  /** alpha is the form alpha[i][j] is a double value for the ith node (indexing from 0) and then the jth parent combination
   * (indexing from 0) for that node **/
  /** the key part of this function is simply to find the correct entry in alpha which has been reqested by g() **/
  
  /** NOTE: we know the node as this is V[dropped] so can go right to the correct row in alpha **/ 
 int i,j;
 int check=0;
 double alpha_max=0.0;
 /*Rprintf("lenV=%d\n",lenV);*/
       /*if(lenV==0){Rprintf("empty\n");}*/
    /* Rprintf("get_alpha_parents()\t");  for(i=0;i<lenV;i++){Rprintf("%d ",V[i]);}Rprintf("\n");*/
       /*for(i=0;i<lenV-1;i++){Rprintf("%d ",V2[i]);}Rprintf("\n");*/
 
 /** IF alpha_i(S) and S has cardinality greater than maxparents then first calc the value for this alpha and then lookup to see where it is - SEE CAVEAT BELOW **/      
 if(lenV-1>maxparents){/** calc the value for this alpha lenV-1 is length of V2 **/
 alpha_max=get_alpha_no_mobius_max(dropped, V2, parents, parents_numparents, nodesid, ptr_score, numRows, numNodes, maxparents, lenV-1);
   for(i=0; i<numRows;i++){/** now find out which row in ptr_score this corresponds to */
     if(ptr_score[i]==alpha_max){/** got a match so return this index as this corresponds to the appropriate row in parents **/
                                 return(i); /** IMPORTANT NOTE: this looks for a value which is an exact match is multiple alphas have the same value then this might
                                                not be correct - although getting an exact match seems highly unlikely **/
    }
    }
   }
 
 
 for(i=0;i<numRows;i++){/** for each possible subset (e.g. alpha) - just need to locate the correct entry in alpha */ 
   /*Rprintf("i=%d\n",i);*/
   /** case 1. empty set **/
   if(   parents_numparents[i]==(lenV-1) /** got right length of subset, -1 as this is analogous to V2 */
      && parents_numparents[i]==0 /** and the empty set **/
      && nodesid[i]==dropped /** got correct node */
   ){ /** we are done **/
      /*Rprintf("alpha[ V[dropped] ][i%(numRows/numNodes)]=%f\n",alpha[ V[dropped] ][i%(numRows/numNodes)]);*/
      return(alphalookup[ dropped ][i%(numRows/numNodes)]);} /** IMPORTANT NOTE - THIS LINE IS NOT THE SAME AS IN get_alpha() as we don't have [ V[dropped] ] etc **/
      
   if(   parents_numparents[i]==(lenV-1) /** -1 as this is analogous to V2 */
      && parents_numparents[i]!=0 
      && nodesid[i]==dropped
    ){ check=0;
       for(j=0;j<lenV-1;j++){/** for each member of V2 **/
	 /*Rprintf("j=%d\n",j);*/
	   
	   if(parents[i][ V2[j] ]==1){check++;}
       }
       /*Rprintf("check=%d\n",check);*/
       if(check==parents_numparents[i]){/** we have both the correct number of parents and the correct combination of parents */
	 /*Rprintf("||");Rprintf("node=%d ",nodesid[i]);for(j=0;j<numNodes;j++){Rprintf("%d ",parents[i][j]);}Rprintf("||\n");
	 Rprintf("alpha[ V[dropped] ][i%(numRows/numNodes)]=%f\n",alpha[ V[dropped] ][i%(numRows/numNodes)]);*/
	 return(alphalookup[ dropped ][i%(numRows/numNodes)]);}
   }
 }
 
 /** if we get to here is it because the number of parents is less than numNodes-1 e.g. a restricted parent set **/
 /** hence some of the alpha_i(S) will not exist and therefore return 0.0 **/
   error("should never get here! - get_alpha_parents()\n");
 return(1);  

  
}
/** *********************************************************************************************/
/** *********************************************************************************************/
/** *********************************************************************************************/
/** *********************************************************************************************/
double get_alpha_no_mobius(int node, int *V2, int **parents, int *parents_numparents, int *nodesid, double *ptr_score, int numRows,int numNodes,int maxparents, int lenV2)
{
 /** want to calculate the sum of all rows in parents where node i = nodesid  and parents are members of the group of subsets of lenV2 with less than maxparents**/
  
  /** V2 contains the sets and its length is lenV-1 **/
 int i,j,k;
 /*  Rprintf("|V2|--");for(i=0;i<lenV2;i++){Rprintf("%d ",V2[i]);}Rprintf("--V2\n");
   Rprintf("node passed=%d\n",node);*/
  /*int *V2tmp=(int *)R_alloc(lenV2,sizeof(int));*/
int  *V2tmp=(int *)malloc(lenV2*sizeof(int));
  int comblength=0;
  int match=0;
  double alpha_value=0.0;
  int maxcardinality;
  gsl_combination *c;
  
 
  /** first step is to determine what the maximum subset size is e.g. it might be less than maxparents **/
  if(lenV2<=maxparents){maxcardinality=lenV2;
    } else {maxcardinality=maxparents;}
 
 /*printf ("All subsets of {0,1,2,3} by size:\n");*/
 
 if(lenV2==0){/** special case - empty set **/ 
              for(j=0;j<numRows;j++){
		 if(  nodesid[j]            == node /** first check nodeid is correct **/
		   && parents_numparents[j] == 0) /** have correct number of parents*/
		 {/*Rprintf("%5.10e\n",ptr_score[j]);*/return(ptr_score[j]);}
   }}
		  
		  
       for (i = 0; i <= maxcardinality; i++)
         {	   
           c = gsl_combination_calloc (lenV2, i);
           do
             {	       
               comblength=gsl_combination_k(c);
               /*for(j=0;j< comblength;j++){Rprintf("%d ",V2[ gsl_combination_data(c)[j] ]);}Rprintf("\n");*/
	       for(j=0;j<lenV2;j++){V2tmp[j]=0;} /** reset V2tmp to empty **/
	       for(j=0;j< comblength;j++){V2tmp[j]=V2[gsl_combination_data(c)[j] ];} /** V2tmp contains the parents */
	       /** have a subset of V2 so sum over all the appropriate rows in ptr_score*/
	       /*Rprintf("want to find this\n");
	       for(k=0;k<comblength;k++){Rprintf("%d ",V2tmp[k]);}Rprintf("\t node=%d \n",node);*/
	       
	       for(j=0;j<numRows;j++){
		 /*for(k=0;k<numNodes;k++){Rprintf("%d ",parents[j][k]);}Rprintf("nodeid=%d numpars=%d\n",nodesid[j],parents_numparents[j]);*/
		 if(  nodesid[j]            == node /** first check nodeid is correct **/
		   && parents_numparents[j] == comblength) /** have correct number of parents*/
		 { match=1;
		   for(k=0;k<comblength;k++){ /** check whether parents match **/
		       if(parents[j][ V2tmp[k] ] == 1){match=1;} else {match=0;break;} 
		                             } 
	           if(match==1){alpha_value+= ptr_score[j];break;}/** got a hit so add this value to the total **/
                }
	       }
	       if(match==0){for(k=0;k<comblength;k++){Rprintf("|%d ",V2tmp[k]);} Rprintf("\t lenV2=%d i=%d\n",lenV2,i);error("no match found - should never get here!");}
             }
           while (gsl_combination_next (c) == GSL_SUCCESS);
           gsl_combination_free (c);
	   
	}
free(V2tmp);	
/*Rprintf("%5.10e\n",alpha_value);*/	
return(alpha_value);	
  
}
/** *********************************************************************************************/
/** *********************************************************************************************/
/** *********************************************************************************************/
/** *********************************************************************************************/
double get_alpha_no_mobius_max(int node, int *V2, int **parents, int *parents_numparents, int *nodesid, double *ptr_score, int numRows,int numNodes,int maxparents, int lenV2)
{
 /** want to calculate the sum of all rows in parents where node i = nodesid  and parents are members of the group of subsets of lenV2 with less than maxparents**/
  
  /** V2 contains the sets and its length is lenV-1 **/
 int i,j,k;
 /*  Rprintf("|V2|--");for(i=0;i<lenV2;i++){Rprintf("%d ",V2[i]);}Rprintf("--V2\n");
   Rprintf("node passed=%d\n",node);*/
 /* int *V2tmp=(int *)R_alloc(lenV2,sizeof(int));*/
 int  *V2tmp=(int *)malloc(lenV2*sizeof(int));
 int comblength=0;
  int match=0;
  double alpha_value;
#ifdef NotLOG 
  alpha_value=0.0;
#endif
#ifdef UseLOG  
  alpha_value= -DBL_MAX;
#endif   
  int maxcardinality;
  gsl_combination *c;
  
 
  /** first step is to determine what the maximum subset size is e.g. it might be less than maxparents **/
  if(lenV2<=maxparents){maxcardinality=lenV2;
    } else {maxcardinality=maxparents;}
 
 /*printf ("All subsets of {0,1,2,3} by size:\n");*/
 
 if(lenV2==0){/** special case - empty set **/ 
              for(j=0;j<numRows;j++){
		 if(  nodesid[j]            == node /** first check nodeid is correct **/
		   && parents_numparents[j] == 0) /** have correct number of parents*/
		 {/*Rprintf("%5.10e\n",ptr_score[j]);*/return(ptr_score[j]);}
   }}
		  
		  
       for (i = 0; i <= maxcardinality; i++)
         {	   
           c = gsl_combination_calloc (lenV2, i);
           do
             {	       
               comblength=gsl_combination_k(c);
               /*for(j=0;j< comblength;j++){Rprintf("%d ",V2[ gsl_combination_data(c)[j] ]);}Rprintf("\n");*/
	       for(j=0;j<lenV2;j++){V2tmp[j]=0;} /** reset V2tmp to empty **/
	       for(j=0;j< comblength;j++){V2tmp[j]=V2[gsl_combination_data(c)[j] ];} /** V2tmp contains the parents */
	       /** have a subset of V2 so sum over all the appropriate rows in ptr_score*/
	       /*Rprintf("want to find this\n");
	       for(k=0;k<comblength;k++){Rprintf("%d ",V2tmp[k]);}Rprintf("\t node=%d \n",node);*/
	       
	       for(j=0;j<numRows;j++){
		 /*for(k=0;k<numNodes;k++){Rprintf("%d ",parents[j][k]);}Rprintf("nodeid=%d numpars=%d\n",nodesid[j],parents_numparents[j]);*/
		 if(  nodesid[j]            == node /** first check nodeid is correct **/
		   && parents_numparents[j] == comblength) /** have correct number of parents*/
		 { match=1;
		   for(k=0;k<comblength;k++){ /** check whether parents match **/
		       if(parents[j][ V2tmp[k] ] == 1){match=1;} else {match=0;break;} 
		                             } 
	           if(match==1){if(ptr_score[j]>alpha_value){alpha_value=ptr_score[j];} /** only update is found a bigger value */
	                        break;}/** got a hit so add this value to the total **/
                }
	       }
	       if(match==0){for(k=0;k<comblength;k++){Rprintf("|%d ",V2tmp[k]);} Rprintf("\t lenV2=%d i=%d\n",lenV2,i);error("no match found - should never get here!");}
             }
           while (gsl_combination_next (c) == GSL_SUCCESS);
           gsl_combination_free (c);
	   
	}
	
free(V2tmp);
/*Rprintf("%5.10e\n",alpha_value);*/	
return(alpha_value);	
  
}



