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
#include "node_binomial.h"
#include "node_gaussian.h"
#include "node_poisson.h"
#include "node_binomial_rv_Rsupp.h"
#include <time.h>

#define DEBUG_12

/** ******************************************************************************************************************************************************/
/** compute all the node scores to create a cache ********************************************************************************************************/                                         
/** ******************************************************************************************************************************************************/
SEXP fitnodes(SEXP R_obsdata, SEXP R_children, SEXP R_nodecache, SEXP R_numVars, SEXP R_numRows, SEXP R_numparspernode, SEXP R_vartype, SEXP R_maxparents,
	      SEXP R_priors_mean, SEXP R_priors_sd,SEXP R_priors_gamshape,SEXP R_priors_gamscale,
			 SEXP R_maxiters, SEXP R_epsabs, SEXP R_verbose, SEXP R_errorverbose, SEXP R_whichnodes,SEXP R_printoutputfreq,
	                 SEXP R_groupedvars, SEXP R_groupids, SEXP R_epsabs_inner,SEXP R_maxiters_inner,
	                 SEXP R_finitestepsize, SEXP R_hparams, SEXP R_maxiters_hessian)
{
  
/** ****************/
/** declarations **/
unsigned int i,j,k,index;
int errverbose,verbose;
datamatrix data,designmatrix;
cache nodecache;
const double priormean=asReal(R_priors_mean);/*Rprintf("priormean=%f %f\n",priormean[0],priormean[5]);*/
const double priorsd=asReal(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double priorgamshape=asReal(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double priorgamscale=asReal(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
int numVars=asInteger(R_numVars);
int numRows=asInteger(R_numRows);
int outfreq=asInteger(R_printoutputfreq);
int numVarsinCache=LENGTH(R_whichnodes);
int *whichnodes=INTEGER(R_whichnodes);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%u ",vartype[i]);} Rprintf("\n");*/
network dag;
/*storage nodescore;*/
SEXP listresults;
SEXP tmplistentry;
int curnode=0;
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
clock_t startclock=0; clock_t endclock=0;double elapsed;int timingon=1;
verbose=asInteger(R_verbose);
errverbose=asInteger(R_errorverbose);

/** end of argument parsing **/

/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
PROTECT(listresults = allocVector(VECSXP, 3));
/*for(i=0;i<3;i++){*//** create TWO vectors to be returned, one for the score, one for
                     an error code for each score **/
  
  PROTECT(tmplistentry=NEW_NUMERIC(numRows));
  SET_VECTOR_ELT(listresults, 0, tmplistentry);
  UNPROTECT(1);
  PROTECT(tmplistentry=NEW_INTEGER(numRows));
  SET_VECTOR_ELT(listresults, 1, tmplistentry);
  UNPROTECT(1);
  PROTECT(tmplistentry=NEW_NUMERIC(numRows));
  SET_VECTOR_ELT(listresults, 2, tmplistentry);
  UNPROTECT(1);
  
				/*}*/
/** *******************************************************************************
***********************************************************************************
 STEP 1. read in combinations of parents into cache - scores all zero'd */

make_nodecache(&nodecache, numVarsinCache, numVars, numRows,R_numparspernode, R_children, R_nodecache, (SEXP)NULL);
/*printCACHE(&nodecache,1);*/

/** create the observed data */
make_data(R_obsdata,&data,R_groupids);
/*printDATA(&data,1);*/

make_dag(&dag, numVars,(SEXP)NULL,1,R_vartype,&maxparents,R_groupedvars);/** create an empty network but with max.parents set **/
/*printDAG(&dag,2);*/

Rprintf("Processing a total of %d node combinations\n",numRows);
/** now neet to iterate over each node in the cache and get its score **/
/*gsl_set_error_handler_off();Rprintf("Warning - turning off GSL Error handler\n");*/ 
 index=0;
 /*for(i=0;i<nodecache.numVars;i++){*//** for each variable */
 for(i=0;i<numVarsinCache;i++){/** for each variable passed */
   
   curnode=whichnodes[i]-1;/** -1 since R indexes start at unity and not zero **/ 
                           /** important NOTE: the nodecache is indexed in terms of i and NOT curnode e.g. it runs from 0,,,numVarsinCache-1 but the actual nodeid are in whichnodes */ 
   for(j=0;j<nodecache.numparcombs[i];j++){/** for each parent combination for this variable */
                         
         /** copy the parent combination in cache into dag  **/
	 for(k=0;k<numVars;k++){dag.defn[curnode][k]=nodecache.defn[i][j][k];}
	   if(timingon && dag.groupedVars[i]){startclock = clock();}
                         switch(dag.varType[curnode])  /** choose which type of node we have */
                         {
			   
			   case 1:{ /** binary/categorical node */
			           if(dag.groupedVars[i]){/** have grouped binary variable so node is a glmm */
				      calc_node_Score_binary_rv_R(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								  0,epsabs_inner,maxiters_inner,finitestepsize, verbose,
								  h_guess,h_epsabs,maxiters_hessian);  
				    } else {/** not grouped so node is a glm **/  
                                    calc_node_Score_binary(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,maxiters,epsabs,0); }
				    REAL(VECTOR_ELT(listresults,0))[index]=dag.nodeScores[curnode]; /** store the score **/
				    INTEGER(VECTOR_ELT(listresults,1))[index]=dag.nodeScoresErrCode[curnode];
				    REAL(VECTOR_ELT(listresults,2))[index]=dag.hessianError[curnode];
				    index++;
                                    break;
                                   }
                         
                           case 2:{ /** gaussian node */
                                    calc_node_Score_gaus(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs, 0);
				    REAL(VECTOR_ELT(listresults,0))[index]=dag.nodeScores[curnode];
				    INTEGER(VECTOR_ELT(listresults,1))[index]=dag.nodeScoresErrCode[curnode];
				    REAL(VECTOR_ELT(listresults,2))[index]=dag.hessianError[curnode]; 
				    index++;
                                    break;
                                   }
			   
			    case 3:{ /** poisson node */
                                    calc_node_Score_pois(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,maxiters,epsabs, 0);
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    REAL(VECTOR_ELT(listresults,0))[index]=dag.nodeScores[curnode];
				    INTEGER(VECTOR_ELT(listresults,1))[index]=dag.nodeScoresErrCode[curnode];
				    REAL(VECTOR_ELT(listresults,2))[index]=dag.hessianError[curnode]; 
                                    index++;
				    break;
                                   }
			        
                           default: {Rprintf("dag.varType[i]=%d\n",dag.varType[curnode]);error("in default switch - should never get here!");}                                          
                         }
                 
                 if(timingon && dag.groupedVars[i]){endclock = clock();                                     
                 elapsed = ((double) (endclock - startclock)) / CLOCKS_PER_SEC;
                 Rprintf("CPU time: %10.6f secs\n",elapsed);
		 }                                     
                         R_CheckUserInterrupt();/** allow an interupt from R console */ 
                         
   if(index%outfreq==0){Rprintf("%d of %d\n",index,numRows);}
   }
 }     

/*gsl_set_error_handler (NULL);*//** restore the error handler*/
free_dag(&dag);/** free any memory not allocated using Ralloc() */
UNPROTECT(1);

return(listresults);

}


