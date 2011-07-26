#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define DEBUG_df_dm_no

/** convert integer data.frame into datamatrix structure for passing to C function */
void df_to_dm(SEXP R_obsdata,datamatrix *obsdata, SEXP R_numVarlevels)
{
int numDataPts,numVars,i,j;
numVars=LENGTH(R_obsdata);/** number of columns in data.frame */
numDataPts=LENGTH(VECTOR_ELT(R_obsdata,0));
int *numVarlevels_ptr=INTEGER_POINTER(R_numVarlevels);
int *numVarlevels, **data, *tmpdata;

/** create a copy of R_vector numVarLevels - could perhap use R internal vector directly but safer this way? */
numVarlevels=(int *)R_alloc( (numVars),sizeof(int));
for(i=0;i<numVars;i++){numVarlevels[i]=numVarlevels_ptr[i];}

/** create a copy of R_data.frame into 2-d C array of ints - note: each CASE is an array NOT each variable */ 
data=(int **)R_alloc( (numDataPts),sizeof(int*));/** number of ROWS*/
	for(i=0;i<numDataPts;i++){tmpdata=(int *)R_alloc( numVars,sizeof(int)); data[i]=tmpdata;} 
 
  for(i=0;i<numDataPts;i++){/** for each CASE/observation **/
     for(j=0;j<numVars;j++){/** for each variable **/
                             data[i][j]=INTEGER(VECTOR_ELT(R_obsdata,j))[i];/** copy data.frame cell entry into C 2-d array entry **/
                              }
     }
 
#ifdef DEBUG_df_fm    
Rprintf("df_to_dm: numVars=%u numDataPts=%u\n",numVars,numDataPts);
for(i=0;i<numVars;i++){Rprintf("numVarLevels=%u\n",numVarlevels[i]);} 

for(i=0;i<numDataPts;i++){for(j=0;j<numVars;j++){Rprintf("%u ",data[i][j]);} Rprintf("\n");}    

#endif

obsdata->data=data;/** original observed data */
obsdata->numVarlevels=numVarlevels;/** number of different categories in each variable */
obsdata->numVars=numVars;/** total number of variables */
obsdata->numDataPts=numDataPts;/** total number of case/observations */

}
/** **********************************************************************************************/
/** **********************************************************************************************/
/** **********************************************************************************/
void store_results(SEXP R_listresults,network *dag, int iter, SEXP ans, int verbose){

int *rans;
int i,j;
/** store the network score **/
REAL(VECTOR_ELT(R_listresults,0))[iter]=dag->networkScore;
/** store the network as a matrix **/
       rans = INTEGER(ans);
       for(i = 0; i < dag->numNodes; i++) {/** fill by row */
         for(j = 0; j < dag->numNodes; j++){
           rans[i + (dag->numNodes)*j] = dag->defn[i][j];}
       }
SET_VECTOR_ELT(R_listresults, iter+1, ans);

if(verbose){
for(i = 0; i < dag->numNodes; i++) {
         for(j = 0; j < dag->numNodes; j++){Rprintf("%d|",dag->defn[i][j]);}Rprintf("\n");}Rprintf("\n");
}

}


