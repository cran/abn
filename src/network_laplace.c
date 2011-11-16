/** this includes supplementary code which is specific to the mixed cts and discrete models using laplace **/
/** laplace still uses some functions in network.c but some needed to be adjusted for mixed of data types **/
#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"        
#include "network_laplace.h"

void build_init_dag_mixed(network *dag,datamatrix *obsdata,int maxparents){
/** build and initialise a matrix which holds the definitions for each node in the network, all nodes are assumed independent initially */
/** use a square form, rows are children and cols parents 
    0 - 00001000110000, node 0 has parents comprising of node 4,8,9
    1 - 00010010011001
    2 - 01010100001000
    etc
*/

int **model_defn,*model_defnA, **model_banlist, *model_banlistA, **model_retainlist, *model_retainlistA, **model_startlist, *model_startlistA;
int i,j;

model_defn=(int **)R_alloc( (obsdata->numVars),sizeof(int*));/** each row is a variable **/ 
	for(i=0;i< obsdata->numVars;i++){model_defnA=(int *)R_alloc( obsdata->numVars,sizeof(int)); model_defn[i]=model_defnA;}

model_banlist=(int **)R_alloc( (obsdata->numVars),sizeof(int*));/** each row is a variable **/ 
	for(i=0;i< obsdata->numVars;i++){model_banlistA=(int *)R_alloc( obsdata->numVars,sizeof(int)); model_banlist[i]=model_banlistA;}

model_retainlist=(int **)R_alloc( (obsdata->numVars),sizeof(int*));/** each row is a variable **/ 
	for(i=0;i< obsdata->numVars;i++){model_retainlistA=(int *)R_alloc( obsdata->numVars,sizeof(int)); model_retainlist[i]=model_retainlistA;}

model_startlist=(int **)R_alloc( (obsdata->numVars),sizeof(int*));/** each row is a variable **/ 
	for(i=0;i< obsdata->numVars;i++){model_startlistA=(int *)R_alloc( obsdata->numVars,sizeof(int)); model_startlist[i]=model_startlistA;}


/** now create an initially null network of independent nodes**/
for(i=0;i<obsdata->numVars;i++){for(j=0;j<obsdata->numVars;j++){model_defn[i][j]=0;}}

/** now create an initially null banlist **/
for(i=0;i<obsdata->numVars;i++){for(j=0;j<obsdata->numVars;j++){model_banlist[i][j]=0;}}

/** send results back to calling environment **/
dag->defn=model_defn;
dag->numNodes=obsdata->numVars;
dag->maxparents=maxparents;
dag->banlist=model_banlist;
dag->retainlist=model_retainlist;
dag->startlist=model_startlist;
}

/** ***************************************************************************/
void init_network_score_mixed(storage *nodescore,network *dag){

int *parentindexes;
parentindexes=(int *)R_alloc(dag->maxparents,sizeof(int));
nodescore->parentindexes=parentindexes;

}

