#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "utility_fns.h"
#include "network.h"
#include <R_ext/Utils.h>
#include "scorereuse.h"
#include <time.h>
#include <gsl/gsl_errno.h>

#define DEBUG_12

SEXP hillsearchfornetwork(SEXP R_obsdata, SEXP R_dag,SEXP R_useK2,SEXP R_maxparents,SEXP R_priorpernode, SEXP R_numVarLevels, 
                      SEXP R_nopermuts, SEXP R_shuffle, SEXP R_labels, SEXP R_nosearches, SEXP R_dag_retain, SEXP R_dag_start,SEXP R_db_size, 
		      SEXP R_localdb, SEXP R_timing, SEXP R_enforce_db_size)
{
/** ****************/
/** declarations **/
unsigned int /*numObs,numNodes,*/i,maxparents,nopermuts,db_size,enforce_db_size;
unsigned int useK2;
double priordatapernode;
datamatrix obsdata;
network dag,dag_scratch,dag_opt1,dag_opt2,dag_opt3,dag_best;
unsigned int nosearches;
unsigned int first;
int iter=0;
unsigned int maxlinkspossible;
struct database prevNodes;/** this will store scores for previous nodes for re-use */
clock_t start=0; clock_t end=0;
double elapsed;
/** end of declarations*/
/** *******************/
/*GetRNGstate();*/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
/*numNodes=LENGTH(R_obsdata);
numObs=LENGTH(VECTOR_ELT(R_obsdata,0));
obsdata.numVars=numNodes;*/
maxparents=asInteger(R_maxparents);
useK2=asInteger(R_useK2);
priordatapernode=asReal(R_priorpernode);
nopermuts=asInteger(R_nopermuts);
nosearches=asInteger(R_nosearches);
db_size=asInteger(R_db_size);
enforce_db_size=asInteger(R_enforce_db_size);
SEXP listresults;
SEXP scorevector;
double lognetworkscore/*,bestlognetworkscore*/;
SEXP ans;
int networkindex=0;
/*int verbose=0;*/
cycle cyclestore;
storage nodescore;
int timingon=asInteger(R_timing);/** whether to turn on timing output or not **/
int localdb=asInteger(R_localdb);/** whether to reset the search database at the start of each search **/
/** end of argument parsing **/

/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/
         PROTECT(listresults = allocVector(VECSXP, nosearches*2+1));/** number of elements in the outer most list 
                                                                        two matrices for each search plus one score vector**/
         PROTECT(scorevector=NEW_NUMERIC(nosearches*2));/** a single vector containing the network score for each step of the search */
         SET_VECTOR_ELT(listresults, 0, scorevector);/** assign the above vector to the first entry in the R list */
         UNPROTECT(1);
/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm(R_obsdata,&obsdata, R_numVarLevels);
/** checked 21/05 - seems to work fine */

/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag(&dag,&obsdata,maxparents);
maxlinkspossible=((dag.numNodes*dag.numNodes)-dag.numNodes);/** used in generate_random_dag()*/
/** all this does is to set internally set dag->banlist[child][parent] etc **/
setbanlist(&dag,R_dag);/** create banned links in initial search graph construction **/
/** NOW for new part - retain list - these arcs must be kept in every new model found */
setretainlist(&dag,R_dag_retain);/** note - sets dag->retainlist*/
/*setstartlist(&dag,R_dag_start);*//** note - sets dag->startlist*/

build_init_dag(&dag_scratch, &obsdata,maxparents);/** create a second dag - a working copy for adding arcs etc **/
setbanlist(&dag_scratch,R_dag);/** create banned links in working copy - **/
setretainlist(&dag_scratch,R_dag_retain);
/*setstartlist(&dag_scratch,R_dag_start);*/

build_init_dag(&dag_opt1, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best add arc etc **/
build_init_dag(&dag_opt2, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best removed arc etc **/
build_init_dag(&dag_opt3, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best reversed arc etc **/

build_init_dag(&dag_best, &obsdata,maxparents); /** simply used to hold the best network found in each search*/

init_network_score(&nodescore,&dag);/** initilise storage for network score **/
init_random_dag(&nodescore,&dag);/** initilise storage for random dag **/
init_hascycle(&cyclestore,&dag); /** initialise storage but needs to be passed down through generate_random_dag etc */
init_nodedatabase(&prevNodes,&dag, db_size,1);/** memory allocation */

for(i=0;i<nosearches;i++){/**out loop for random re-start hill climber */
if(timingon){start = clock();}  
if(localdb){init_nodedatabase(&prevNodes,&dag, db_size,0);}/** resets the database but leave memory alone **/
Rprintf("search number...%d\n",i);
generate_random_dag(&cyclestore,&nodescore,&dag,nopermuts,maxparents,R_shuffle,i*maxlinkspossible,i, R_dag_start); /** replace the dag->defn with a new structure **/ 
calc_network_Score_reuse(&nodescore,&dag,&obsdata,priordatapernode, useK2,0,R_labels,&prevNodes,enforce_db_size);/** 0 is to turn off printing out parameters for each node */

Rprintf("initial network: (log) network score = %f\n",dag.networkScore);

copynetworkdefn(&dag,&dag_best);/** make a copy of the initial network in case new network is worse*/
dag_best.networkScore=dag.networkScore;/** make a copy of the initial score in case new score is worse*/

/** creates an R matrix which will contain the network structure of the initial random network 
             store_results() sets this matrix into the outer R list and also sets the network score for this in the vector in the list*/
         PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
         store_results(listresults,&dag_best,networkindex++, ans,0);/** first 0 is the index of networks to save ignoring first scorevector */
	 UNPROTECT(1);
         /** we have now stored the structure and network score of the initial network */

/** now for an individual stepwise search given the initial random network*/
lognetworkscore=dag.networkScore;/** start off with score of the random starting network */
           first=1;iter=1;
           while(dag.networkScore>lognetworkscore || first){
                lognetworkscore=dag.networkScore; 
                /** do next step in iterative single arc add/remove/reversal search */
                hill_climb_iter(&nodescore,&cyclestore,&dag, 
                                &dag_scratch,
                                &dag_opt1,
                                &dag_opt2,
                                &dag_opt3,
                                maxparents,&obsdata, priordatapernode,useK2,0,R_labels,&prevNodes,enforce_db_size);/** &dag will have new best network*/
                R_CheckUserInterrupt();/** allow an interupt from R console */ 
                /** got a better network then update score and structure, if not do nothing and while() will terminate */
                if(dag.networkScore>lognetworkscore){
                  copynetworkdefn(&dag,&dag_best);/** copy new network */
                  dag_best.networkScore=dag.networkScore;/** copy new network's score */
                  
                  /*Rprintf("search iteration...%d new score=%f\n",iter,dag.networkScore);*/
                  iter++;
                  }
                first=0;/** flag used on first iteration - now unset */
                }
Rprintf("best network: (log) network score = %f\n",dag_best.networkScore); 
          /** now got the best network from current search so save back to R **/         
          PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
                           store_results(listresults,&dag_best,networkindex++,ans,0);
          UNPROTECT(1);                 
          
if(timingon){end = clock();
 elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
 Rprintf("CPU time:%10.6f\n",elapsed);}
 
} /** END OF outer for loop **/

/** some diagnostic messages */
 if(prevNodes.nodecacheexceeded){Rprintf("Note: db.size of %u exceeded. Increasing this may be more efficient\n",prevNodes.length);}        

/*UNPROTECT(nosearches*2+1+1); */
UNPROTECT(1);

/*free_dag(&dag);
free_dag(&dag_scratch);
free_dag(&dag_opt1);
free_dag(&dag_opt2);
free_dag(&dag_opt3);
free_dag(&dag_best);
*/
return(listresults);

}



