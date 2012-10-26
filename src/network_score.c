#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h" 
#include "network_score.h"
#include "node_binomial.h"
#include "node_binomial_rv_Rsupp.h"
#include "node_gaussian.h"
#include "node_poisson.h"
#include "node_binomial_marginals_rv_Rsupp.h"
#include <time.h>
/** ***************************************************************************************************/
/**  pass a DAG and call the appropriate node score function depending on the distibution of the node */
/** ***************************************************************************************************/
void calc_network_Score(network *dag,datamatrix *obsdata, datamatrix *designmatrix,
				const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs, int verbose, const int errverbose, SEXP results, int storeModes, 
			        double epsabs_inner, int maxiters_inner, double finitestepsize,
			         double h_guess, double h_epsabs, int maxiters_hessian)
{
 int i;
 int numnodes=dag->numNodes;
 clock_t startclock=0; clock_t endclock=0;double elapsed;int timingon=1;
 
 for(i=0;i<numnodes;i++){Rprintf("processing node %i\n",i+1);
 /*i=3;*/
                         if(timingon && dag->groupedVars[i]){startclock = clock();}
                         switch(dag->varType[i])  /** choose which type of node we have */
                         {
                           case 1:{ /** binary/categorical node */ 
			            if(dag->groupedVars[i]){/** have grouped binary variable so node is a glmm */
				      calc_node_Score_binary_rv_R(dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								  storeModes,epsabs_inner,maxiters_inner,finitestepsize,verbose,
								  h_guess,h_epsabs, maxiters_hessian);  
				    } else {/** not grouped so node is a glm **/  
                                      calc_node_Score_binary(dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,maxiters,epsabs,storeModes); 
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    }
				    storeResults(results,dag,storeModes,i,dag->varType[i]);
                                    /*if(verbose){Rprintf("Binary node=%d score=%f\n", i,dag->nodeScores[i]);} */                                                                                        
                                    break;
                                   }
                         
                           case 2:{ /** gaussian node */ 
			           if(dag->groupedVars[i]){/** have grouped gaussian variable so node is a glmm */
				      error("Gaussan node with random effects not yet implemented!\n");
				    } else {/** not grouped so node is a glm **/  
                                      calc_node_Score_gaus(dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs, storeModes);
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    }
				    storeResults(results,dag,storeModes,i,dag->varType[i]);				 
                                    /*if(verbose){Rprintf("Gaussian node=%d score=%f\n",i,REAL(VECTOR_ELT(results,0))[i]);}*/
                                    break;
                                   }
                           
			   case 3:{ /** poisson node */
			           if(dag->groupedVars[i]){/** have grouped poisson variable so node is a glmm */
				      error("Poisson node with random effects not yet implemented!\n");
				    } else {/** not grouped so node is a glm **/  
                                      calc_node_Score_pois(dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,maxiters,epsabs, storeModes);
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    }
				    storeResults(results,dag,storeModes,i,dag->varType[i]);				 
                                    /*if(verbose){Rprintf("Poisson node=%d score=%f\n",i,REAL(VECTOR_ELT(results,0))[i]);}*/
                                    break;
                                   }        
                                   
                                   
                           default: {error("in default switch - should never get here!");}                                          
                         }
                 
                 if(timingon && dag->groupedVars[i]){endclock = clock();                                     
                 elapsed = ((double) (endclock - startclock)) / CLOCKS_PER_SEC;
                 Rprintf("CPU time: %10.6f secs\n",elapsed);
		 }
		 R_CheckUserInterrupt();/** allow an interupt from R console */ 
                         /*if(verbose){Rprintf("individual node score=%f\n",indnodescore);}*/
                        /* lognetworkscore+=indnodescore;*/
}      
                    
 /*dag->networkScore=lognetworkscore;*/
 
            
 }
/** *************************************************************************************************/
/** *************************************************************************************************/
/** ***************************************************************************************************/
/**  pass a DAG and call the appropriate node score function depending on the distibution of the node */
/** ***************************************************************************************************/
void calc_parameter_marginal(network *dag,datamatrix *obsdata, datamatrix *designmatrix,
				const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs, int verbose, const int errverbose, 
			      double *denom_modes, int childid, int paramid,  
			     double epsabs_inner, int maxiters_inner, double finitestepsize,
			     double h_guess, double h_epsabs,int maxiters_hessian,
			     double betafixed, double mlik, double *posterior)
{
  
                         switch(dag->varType[childid])  /** choose which type of node we have */
                         {
                           case 1:{ /** binary/categorical node */
			            if(dag->groupedVars[childid]){/** have grouped binary variable so node is a glmm */
				      calc_binary_marginal_rv_R(dag,obsdata,childid,errverbose, designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								epsabs_inner,maxiters_inner,finitestepsize,verbose,
								h_guess,h_epsabs,maxiters_hessian,
								denom_modes, paramid, betafixed,mlik, posterior); 
			       
				    } else {/** not grouped so node is a glm **/  
                                      calc_binary_marginal(dag,obsdata,childid,errverbose, designmatrix, priormean, priorsd,maxiters,epsabs,denom_modes,paramid,betafixed, mlik, posterior
							 ); 
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    }
				    
                                    
				    
                                    /*if(verbose){Rprintf("Binary node=%d score=%f\n", i,REAL(results)[i]);}*/                                                                                        
                                    break;
                                   }
                        
                           case 2:{ /** gaussian node */
			            
                                    calc_gaussian_marginal(dag,obsdata,childid,errverbose, designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,
						       epsabs,denom_modes,paramid,betafixed, mlik, posterior);
				    			 
                                    /*if(verbose){Rprintf("Gaussian node=%d score=%f\n",i,REAL(results)[i]);}*/
                                    break;
                                   }
                                   
                                   
                            case 3:{ /** poisson node */
			            
                                    calc_poisson_marginal(dag,obsdata,childid,errverbose, designmatrix, priormean, priorsd,maxiters,epsabs,denom_modes,paramid,betafixed, mlik, posterior
							 ); 
				   
                                    /*if(verbose){Rprintf("Binary node=%d score=%f\n", i,REAL(results)[i]);}*/                                                                                        
                                    break;
                                   }
                                   
                           default: {error("in default switch - should never get here!");}                                          
                         }                                     
                         R_CheckUserInterrupt();/** allow an interupt from R console */ 
                         /*if(verbose){Rprintf("individual node score=%f\n",indnodescore);}*/
                        /* lognetworkscore+=indnodescore;*/
 /*}*/      /*}*/ 
                    
 /*dag->networkScore=lognetworkscore;*/
 
            
 }
/** *************************************************************************************************/
/** *************************************************************************************************/
/** copy results from C into R SEXP *****************************************************************/ 
void storeResults(SEXP results,network *dag,int storeModes, int nodeid, int vartype)
{
  int i;
/** format is results has one vector for each node of length numNodes+3, the format is 
    vec= node score, then each of the parameter modes - including entries DBL_MAX if node not in model**/  
  REAL(VECTOR_ELT(results,nodeid))[0]=dag->nodeScores[nodeid]; /** first entry is the mlik */
  REAL(VECTOR_ELT(results,nodeid))[1]=dag->nodeScoresErrCode[nodeid];
  REAL(VECTOR_ELT(results,nodeid))[2]=dag->hessianError[nodeid];
  /** next get all the parameters including intercept and extra term for gaussian precision parameter etc**/
  /** +4 = mlik+intercept+gauss precision + group precision**/
  for(i=1+2;i<dag->numNodes+2+4;i++){REAL(VECTOR_ELT(results,nodeid))[i]=gsl_matrix_get(dag->modes,nodeid,i-1-2);}
  
  
}

