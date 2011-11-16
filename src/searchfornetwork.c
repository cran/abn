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


#define DEBUG_12

SEXP searchfornetwork(SEXP R_obsdata, SEXP R_dag,SEXP R_useK2,SEXP R_maxparents,SEXP R_priorpernode, SEXP R_numVarLevels, 
                      SEXP R_nopermuts, SEXP R_shuffle, SEXP R_labels, SEXP R_dag_retain, SEXP R_dag_start, SEXP R_db_size, SEXP R_enforce_db_size)
{
/** ****************/
/** declarations **/
unsigned int i,maxparents,nopermuts,db_size,enforce_db_size;/*,numObs,numNodes*/;
unsigned int useK2;
double priordatapernode;
datamatrix obsdata;
network dag,dag_scratch,dag_opt1,dag_opt2,dag_opt3,dag_best;
/*unsigned int nosearches=1;*//** HARD CODED SINCE ONLY SINGLE SEARCH CONDUCTED **/
unsigned int first;
int iter=0;
struct database prevNodes;/** this will store scores for previous nodes for re-use */
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
db_size=asInteger(R_db_size);
enforce_db_size=asInteger(R_enforce_db_size);
SEXP listresults=0;/* just to avoid uninitialised erorr */
SEXP tmplistentry;
double lognetworkscore;/*,bestlognetworkscore;*/
SEXP ans;
int listsize=0;
/*int verbose=0;*/
cycle cyclestore;
storage nodescore;
/** end of argument parsing **/

/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double 
   - currently overkill but useful template **/

/** *******************************************************************************
***********************************************************************************
 STEP 1. convert data.frame in R into C data structure for us with BGM functions */

/** convert integer data.frame into datamatrix structure for passing to C function */
df_to_dm(R_obsdata,&obsdata, R_numVarLevels);
/** checked 21/05 - seems to work fine */

/** initalise network structure - storage for network definition and all (hyper)parameters, 
     this covers any valid network with <=max.parents **/
build_init_dag(&dag,&obsdata,maxparents);
/** all this does is to set internally set dag->banlist[child][parent] etc **/
setbanlist(&dag,R_dag);/** create banned links in initial search graph construction - default is not banned arcs**/
/** NOW for new part - retain list - these arcs must be kept in every new model found */
setretainlist(&dag,R_dag_retain);/** note - sets dag->retainlist*/
/*setstartlist(&dag,R_dag_start);*//** note - sets dag->startlist*/

build_init_dag(&dag_scratch, &obsdata,maxparents);/** create a second dag - a working copy for adding arcs etc **/
setbanlist(&dag_scratch,R_dag);/** create banned links in working copy - again default is no banned arcs**/
setretainlist(&dag_scratch,R_dag_retain);
/*setstartlist(&dag_scratch,R_dag_start);*/

build_init_dag(&dag_opt1, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best add arc etc **/
build_init_dag(&dag_opt2, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best removed arc etc **/
build_init_dag(&dag_opt3, &obsdata,maxparents);/** create a working copy local to hill_climb_iter for holding best reversed arc etc **/

build_init_dag(&dag_best, &obsdata,maxparents); /** simply used to hold the best network found */

init_network_score(&nodescore,&dag);/** initilise storage for network score **/
init_random_dag(&nodescore,&dag);/** initilise storage for random dag **/
init_hascycle(&cyclestore,&dag); /** initialise storage but needs to be passed down through generate_random_dag etc */
init_nodedatabase(&prevNodes,&dag,db_size,1);/** memory allocation */

for(i=0;i<2;i++){/** This is used to run the search TWICE - inefficient but easiest way to deal with R memory allocation
                     since the number of steps taken in the stepwise search are not known in advance and so the list length
                     back to R is of variable length - simple fix is to run the search twice, identical each time, so when i==0
                     this search does all the numerics but does no R storage, we then know how much storage is needed and 
                     repeat search this time starting off the making the necessary R allocations. Inefficient but search a single
                     network does not take long anyway. The if(i==1) are just the storage bits*/

/*if(i==1){verbose=1;}*/ /** this just turns on some output to STDOUT in hill_climb_iter **/

/** THIS PART NEED TO ADJUST *******/
generate_random_dag(&cyclestore,&nodescore,&dag,nopermuts,maxparents,R_shuffle,0,0, R_dag_start); /** replace the dag->defn with a new structure **/ 
calc_network_Score_reuse(&nodescore,&dag,&obsdata,priordatapernode, useK2,0,R_labels,&prevNodes,enforce_db_size);

if(i==1){Rprintf("initial network: (log) network score = %f\n",dag.networkScore);
         /** now arrange the R memory structures **/
         PROTECT(listresults = allocVector(VECSXP, listsize));/** number of elements in the outer most list **/
         PROTECT(tmplistentry=NEW_NUMERIC(listsize-1));/** a single vector containing the network score for each step of the search */
         SET_VECTOR_ELT(listresults, 0, tmplistentry);/** assign the above vector to the first entry in the R list */
				 UNPROTECT(1);
         }

copynetworkdefn(&dag,&dag_best);/** make a copy of the current network in case new network is worse*/
dag_best.networkScore=dag.networkScore;/** make a copy of the current score in case new score is worse*/

if(i==1){/** creates an R matrix which will contain the network structure of the initial random network 
             store_results() sets this matrix into the outer R list and also sets the network score for this in the vector in the list*/
         PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
         store_results(listresults,&dag_best,0, ans,0);
         UNPROTECT(1);
         /** we have now stored the structure and network score of the initial network */
         }

/** now for the actual search process */
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
                                maxparents,&obsdata, priordatapernode,useK2,0,R_labels, &prevNodes,enforce_db_size);/** &dag will have new best network*/
                R_CheckUserInterrupt();/** allow an interupt from R console */ 
                /** got a better network then update score and structure, if not do nothing and while() will terminate */
                if(dag.networkScore>lognetworkscore){
                  copynetworkdefn(&dag,&dag_best);/** copy new network */
                  dag_best.networkScore=dag.networkScore;/** copy new network's score */
                  if(i==1){/** as with initial random network store the structure and the score */
                           PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
                           store_results(listresults,&dag_best,iter,ans,0);          
                           UNPROTECT(1);
                  Rprintf("search iteration...%d new score=%f\n",iter,dag.networkScore);}
                  iter++;
                  }
                first=0;/** flag used on first iteration - now unset */
                }
          listsize=iter+1;/** used to set length of R outer list - checked and this works so dont change */
          
          
   } /** END OF outer for loop which runs the search twice **/
 
/** some diagnostic messages */
 if(prevNodes.nodecacheexceeded){Rprintf("Note: db.size of %u exceeded. %u calls to db after limit reached (n.b. the same nodes may be called multiple times)\n",prevNodes.length,prevNodes.overflownumentries);}  
 
/*UNPROTECT(listsize+1);*/
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



