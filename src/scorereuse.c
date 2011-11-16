#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"        
#include "network.h"
#include <search.h>

#define nolibcSearch
/** IMPORTANT note - using nolibcSearch appears to be quite a bit faster than using the lfind() libc function **/

#if defined nolibcSearch
/** *************************************************************************************/
/** setup memory storage for storing previously calculated nodes **/
/** *************************************************************************************/
void init_nodedatabase(struct database *prevNodes, network *dag, int db_size, int allocMem)
{
 int **knownnodes=0; int *knownnodesA=0; int i;
 double *nodescores=0;
  if(allocMem){/** if allocMem is false then this function just resets the counters and the db is overwritten */
  nodescores=(double *)R_alloc(db_size,sizeof(double));    
  knownnodes=(int **)R_alloc(db_size,sizeof(int*));/** each row is a node**/ 
	for(i=0;i< db_size;i++){knownnodesA=(int *)R_alloc(dag->numNodes+1,sizeof(int)); knownnodes[i]=knownnodesA;} /** +1 as first entry denote the child !!**/  
  
  prevNodes->knownnodes=knownnodes; /** this will hold the array describing the nodes previously found */
  prevNodes->knownscores=nodescores; /** this will hold the scores previously found */
  prevNodes->length=db_size;
  }
  prevNodes->numentries=0;/*** this will keep a count of the number of entries actually used */
  prevNodes->nodecacheexceeded=0;/** a flag - true 1 - if current cache is exceeded */
  prevNodes->overflownumentries=0; /** will hold hold how many times the cache was called after it was full **/
}

/** *************************************************************************************************/
/** searches the "database" for previously calculated nodes */
/** *************************************************************************************************/
int nodescoreisknown(network *dag,int nodeid, double *nodescore, struct database *prevNodes)
{
 unsigned int i,j;              
 unsigned int num_entries=prevNodes->numentries;
 int **knownnodes=prevNodes->knownnodes;
 double *knownscores=prevNodes->knownscores;
 int match=1;
 
 /** step 1. find out whether we already have a score for this node */ 
 /*Rprintf("passed node %u\t",nodeid);
 for(j=0;j<dag->numNodes;j++){Rprintf("%u",dag->defn[nodeid][j]);}Rprintf("\n");
 */ 

 /*lfind (const void *key, void *base, size_t *nmemb, size_t size, comparison_fn_t compar)*/
 /*Rprintf("num_entries=%d\n",num_entries);*/
 for(i=0;i<num_entries;i++){/** each each "currrent" entry **/
      match=1;
      /*for(j=0;j<dag->numNodes;j++){Rprintf("%u",knownnodes[i][j]);}Rprintf("\n"); */
      
      if(knownnodes[i][0]!=nodeid){continue;}/** comparing with correct child/node if not matched then skip to end of out for() **/
      for(j=0;j<dag->numNodes;j++){/** for each node in the current entry check whether its the same as the one passed */
                 if(dag->defn[nodeid][j]!=knownnodes[i][j+1]){match=0; /*found a discrepancy so stop searching this node and move to next one */
                                                       break;}         /** the j+1 is because the first entry in each row is the node id */
                 }
      
      if(match==1){/** then we have found the node we need so set the score and finish**/           
                   /*Rprintf("got existing node \n");*/
                   /*Rprintf("known\n");*/
                   *nodescore=knownscores[i];
                   /*Rprintf("scree=%f %u\n",knownscores[i],i); */ 
                   return(1); }
      }
               
  /** if we get to here we never found any matches and so return false */    
                 
return 0;

}
/** *************************************************************************************************/
/** searches the "database" for previously calculated nodes */
/** *************************************************************************************************/
void storenodescore(network *dag,int nodeid, double nodescore, struct database *prevNodes, int enforce_db_size)
{
 unsigned int j;              
 unsigned int cur_length=prevNodes->length;
 unsigned int num_entries=prevNodes->numentries;
 int **knownnodes=prevNodes->knownnodes;
 double *knownscores=prevNodes->knownscores;
 
 if(prevNodes->numentries<cur_length){/** then do not need to add more memory before adding new record to score **/
                                       
 /** append/copy node onto end of existing num_entires **/
 knownnodes[num_entries][0]=nodeid;/** set the node id */
 for(j=0;j<dag->numNodes;j++){knownnodes[num_entries][j+1]=dag->defn[nodeid][j];} /** j+1 +1 offser since first entry is row id */
                              knownscores[num_entries]=nodescore;
                              /*Rprintf("adding...%u\n",num_entries);*/
 prevNodes->numentries=num_entries+1;/** increment row counter **/
 } else {if(enforce_db_size){error("db.size of %d exceeded. Terminating search. Increase db.size or turn off enforce.db.size\n",cur_length);}
         prevNodes->nodecacheexceeded=1;/** a flag to denote if cache should be increased */
         /** increment row counter here so we know how big to make it next time*/
         prevNodes->overflownumentries++; 
         }
         /** have not yet written code to dynamically increase array as this may be quite cumbersome and inefficient...? */
        

}
#else
/** new search stuff is down here *****************************************************************/
/** *************************************************************************************/
/** setup memory storage for storing previously calculated nodes **/
/** *************************************************************************************/
void init_nodedatabase(struct database *prevNodes, network *dag, int db_size, int allocMem)
{
 int **knownnodes=0; int *knownnodesA=0; int i;
 double *nodescores=0;
 if(allocMem){/** if allocMem is false then this function just resets the counters and the db is overwritten */
  nodescores=(double *)R_alloc(db_size,sizeof(double));    
  knownnodes=(int **)R_alloc(db_size,sizeof(int*));/** each row is a node**/ 
	for(i=0;i< db_size;i++){knownnodesA=(int *)R_alloc(dag->numNodes+2,sizeof(int)); knownnodes[i]=knownnodesA;} 
	/** +2 as first entry denotes an index - row number - and the second is the child node - needed in the lfind routine **/  
 	   
  prevNodes->knownnodes=knownnodes; /** this will hold the array describing the nodes previously found */
  prevNodes->knownscores=nodescores; /** this will hold the scores previously found */
  prevNodes->length=db_size;
  }
  prevNodes->numentries=0;/*** this will keep a count of the number of entries actually used */
  prevNodes->nodecacheexceeded=0;/** will hold how many times the current cache is exceeded */
}
/** *************************************************************************************/
/** *********** an array search comparison function ************************/
/** * this has been checked 04-10-2011, it is not as trivial as it seems due to
the 2-d array of pointers to pointers etc and was devised with much messing about*/
/** current written to compare rows of a 2-d array ignoring the first two entries in each row 
 *  a match is if the two rows are identical */
/** *************************************************************************************/
/** the below struct is needed in nodescoreisknown in order to make the lfind work **/
struct fnetwork {
      int *defn;/** each row a variable and each col the parents of the variable, indexes from 0*/
      int nodeid;
      unsigned int len;
};


int compare_array (const void *key, const void *dbdata)
     { 
      
       int i;
       struct fnetwork *mykey=(struct fnetwork *)key;
       int *v1=mykey->defn; 
       int nodeid=mykey->nodeid;
       int len=mykey->len;/** number of nodes*/
       const int **mydbdata=(const int**)dbdata;
       const int *v2=*mydbdata;/** this is a bit odd - cast on cast - but seems to work ? */
       int curnodeid;
      
      if(v2==0){return(0);} /** this IS needed to avoid seg fault if not match found **/
      curnodeid=v2[1];/** the child node in this row */
      /** some useful debugging printing stuff **/ 
      /*for(i=0;i<len+2;i++){Rprintf("|%u %u|",*(*(mydbdata)+i), v2[i]);}Rprintf("\n");
      Rprintf("testdata=\n");
      for(i=0;i<len;i++){Rprintf("|%u|",v1[i]);}Rprintf("\n");
      */
     if(nodeid!=curnodeid){return 1; /** different child node so no match**/
      }

      /** this is the key part - is two array entries differ then return false n.b. =1 here */ 

      for(i=0;i<len;i++){ /** we ignore the first entry in the comparison **/
             if(v1[i]!=v2[i+2]){/* offset is due to first two entries are an index and child node*/
                               return(1);            }} /* got a discord so finish */
       /** if get to here then must have a match */
      return 0; /** n.b. 0 is a MATCH */
       
     }
/** *********** end of function ************************/
/** *************************************************************************************************/
/** searches the "database" for previously calculated nodes */
/** *************************************************************************************************/
int nodescoreisknown(network *dag,int nodeid, double *nodescore, struct database *prevNodes)
{             
size_t num_entries=prevNodes->numentries;
 int **knownnodes=prevNodes->knownnodes;
 double *knownscores=prevNodes->knownscores;
 int **livehere;
 struct fnetwork mynet;
 
 if(num_entries==0){return(0);}
 /** step 1. find out whether we already have a score for this node */ 
 mynet.defn=dag->defn[nodeid]; /**store pointer to current node definition */ 
 mynet.nodeid=nodeid; /** also store child node id */
 mynet.len=dag->numNodes;/** number of nodes in network **/
 
 /*Rprintf("num_entries=%d nodeid=%d\n",num_entries, nodeid);*/
 
 livehere=(int **)lfind ( &mynet, knownnodes, (size_t*)&num_entries, sizeof(int*), compare_array);

 if(livehere){/** have a hit **/
                   /* printf("got hit - first entry in array is ==>%d\n",*(*livehere+1));*//** return 0th entry in hit **/
                    /** now return correct score **/
		      *nodescore= knownscores[ *(*livehere+0) ];/** [] is the first entry in the array which was a hit **/
        return(1);
 
 } else {/** no hit **//*printf("no hit\n");*/
         return(0);
   }

}
/** *************************************************************************************************/
/** searches the "database" for previously calculated nodes */
/** *************************************************************************************************/
void storenodescore(network *dag,int nodeid, double nodescore, struct database *prevNodes)
{
 unsigned int j;              
 unsigned int cur_length=prevNodes->length;
 size_t num_entries=prevNodes->numentries;
 int **knownnodes=prevNodes->knownnodes;
 double *knownscores=prevNodes->knownscores;
 
 if(prevNodes->numentries<cur_length){/** current number of entries is less than number of slots available **/
                                       
 /** append/copy node onto end of existing num_entires **/
 knownnodes[num_entries][0]=num_entries;/** increment the "row" number */
 knownnodes[num_entries][1]=nodeid;/** set the node id */
 for(j=0;j<dag->numNodes;j++){knownnodes[num_entries][j+2]=dag->defn[nodeid][j];} /** j+1 +1 offser since first entry is row id */
                              knownscores[num_entries]=nodescore;
                              /*Rprintf("adding...%u\n",num_entries);*/
			      
 prevNodes->numentries=num_entries+1;/** increment row counter **/
 } else {prevNodes->nodecacheexceeded=1;/** a flag to denote if cache should be increased */
         }
         /** have not yet written code to dynamically increase array as this may be quite cumbersome and inefficient...? */
        

}
#endif
