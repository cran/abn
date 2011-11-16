#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"        
#include "network.h"
#include "scorereuse.h"

#define DEBUG_dagdefn_no
#define DEBUG_LEARN_NODE_no
#define DEBUG_dagbanlist_no
#define STRINGLENGTH 1000


void build_init_dag(network *dag,datamatrix *obsdata,int maxparents){
/** build and initialise a matrix which holds the definitions for each node in the network, all nodes are assumed independent initially */
/** use a square form, rows are children and cols parents 
    0 - 00001000110000, node 0 has parents comprising of node 4,8,9
    1 - 00010010011001
    2 - 01010100001000
    etc
*/

int **model_defn,*model_defnA, **model_banlist, *model_banlistA,**model_retainlist, *model_retainlistA,**model_startlist, *model_startlistA;
int i,j;
int *numvarlevels=obsdata->numVarlevels;/** hold the number categories for each variable in the data **/
unsigned int *topKLevels=(unsigned int*) R_alloc (maxparents,sizeof(unsigned int));/** will hold the number of levels in the top maxparents number of variables **/ 
int totalMaxParentCombinations=1;
double **nodesparameters,*dirichletparms;
int **nodesparameters_lookup,*parentcomb;
int *nodeparentcombs;

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

/** finally create a vector which will store the number of categories in each parentset per node **/
nodeparentcombs=(int *)R_alloc( obsdata->numVars,sizeof(int)); 

/** send results back to calling environment **/
dag->defn=model_defn;
dag->numNodes=obsdata->numVars;
/*dag->namesNodes=obsdata->namesVars;*//** note this is simply a copy of the address not a copy of data */
dag->numNodeLevels=nodeparentcombs;/** array of the total number of combinations of parent which each node has **/
dag->maxparents=maxparents;
dag->banlist=model_banlist;
dag->retainlist=model_retainlist;
dag->startlist=model_startlist;
/** have got the network definition so now create storage for model parameters at each potential node*/
/** create a 2d array with X rows where X is the max. number of parent combinations possible given
    maxparents and the number of levels each variable could have **/

/** NOTE using nodeparentcombs as temporary array since this has nothing in it yet **/
for(i=0;i<obsdata->numVars;i++){nodeparentcombs[i]=numvarlevels[i];}/** copy into here as R_isort over-rights original array **/
R_isort(nodeparentcombs,obsdata->numVars);/** sort ascending */                                        
for(i=(obsdata->numVars-1);i>=(obsdata->numVars-maxparents)+1-1;i--){
                         topKLevels[(obsdata->numVars-1)-i]=nodeparentcombs[i];} /** NOTE: topKLevels index must go from 0 upwards! */

/** re-set nodeparentcombs **/
for(i=0;i<obsdata->numVars;i++){nodeparentcombs[i]=0;}

for(i=0;i<maxparents;i++){
                          totalMaxParentCombinations*=topKLevels[i];
                          }

dag->maxParentCombinations=totalMaxParentCombinations;

/** now create a 3d structure which will hold all the parameters for a node in the network, big enough for any node **/
/** note that as we only ever learn a model node by node we only need to create storage for one node not all the network 
    nodeparameters[j][k]=parameters for parent combination j, k is the kth dirichlet parameter for this parent combination 
    nodeparameters_lookup[j][k]=parent combination j, k is the category level of one of the nodes parents

    we assume that each node has the maxnumber of parents and then parents with the most number of levels,
    we also assume that the node itself has the most number of possible levels so enough space is reserved
    for any eventuality

    We also create a similar structure but with the final dimension that of the actual parent combination
**/

nodesparameters=(double **)R_alloc( totalMaxParentCombinations,sizeof(double*)); 
                                      for(j=0;j<totalMaxParentCombinations;j++){
				      		dirichletparms=(double *)R_alloc( topKLevels[0],sizeof(double));
                                      		nodesparameters[j]=dirichletparms;	
										}
                                      
/** as above but now the final dimension is a lookup of the actual combination of parent levels corresponding to the parent comb index*/

nodesparameters_lookup=(int **)R_alloc( totalMaxParentCombinations,sizeof(int*)); 
                                      for(j=0;j<totalMaxParentCombinations;j++){
				      		parentcomb=(int *)R_alloc( maxparents,sizeof(int));
                                      		nodesparameters_lookup[j]=parentcomb;}


/** Hence nodeparameters[j][k] = is the kth dirichlet parameter in P(A|B,C) ~ dirichlet(p1,p2,...pk,..)
          nodeparameters_lookup[j][k] = is the observed level of the kth parent variable in P(A|x1,x2,...xj,..) */ 

/** now pass the address back to calling environment **/


dag->nodesparameters=nodesparameters;
dag->nodesparameters_lookup=nodesparameters_lookup;


}

/** *************************************************************************************/
/** read in a DAG definition from an R matrix and copy into C array in network struct **/
/** *************************************************************************************/
void get_dag(network *dag,SEXP R_dag){

unsigned int i,j;

/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
/** dag->defn memory allocated in build_init_dag() */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->defn[i][j]=INTEGER(R_dag)[i+j*dag->numNodes];}} 

#ifdef DEBUG_dagdefn   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->defn[i][j]);} Rprintf("\n");}    
#endif

}
/** *************************************************************************************/
/** read in a DAG definition from a list of R matrices and copy into C array in network struct **/
/** *************************************************************************************/
void get_dag_list(network *dag,SEXP R_dag_list,unsigned int index){

unsigned int i,j;
/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
/** dag->defn memory allocated in build_init_dag() */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->defn[i][j]=INTEGER(VECTOR_ELT(R_dag_list,index))[i+j*dag->numNodes];}} 

#ifdef DEBUG_dagdefn   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->defn[i][j]);} Rprintf("\n");}    
#endif

}

/** **************************************************************************************************/
/** create storage for doing cycle check - to replace previous static allocation calls */
/** **************************************************************************************************/
void init_hascycle(cycle *cyclestore,network *dag){

unsigned int i;
unsigned int numnodes=dag->numNodes;
unsigned int *isactive, *incomingedges;
unsigned int **graph,*graphtmp;
          
   isactive=(unsigned int *)R_alloc(numnodes,sizeof(unsigned  int));
   incomingedges=(unsigned  int *)R_alloc(numnodes,sizeof(unsigned  int));
   graph=(unsigned  int **)R_alloc(numnodes,sizeof(unsigned  int*));/** create storage for a copy of the dag->defn[][] */
   for(i=0;i<numnodes;i++){graphtmp=(unsigned  int *)R_alloc( numnodes,sizeof(unsigned  int)); graph[i]=graphtmp;}

cyclestore->isactive=isactive;
cyclestore->incomingedges=incomingedges;
cyclestore->graph=graph;


}

/** *************************************************************************************************/
/** check for cycle in graph  - do this by checking for a topolgical ordering */
/** *************************************************************************************************/
unsigned int hascycle(cycle *cyclestore,network *dag){

unsigned int numnodes=dag->numNodes;
unsigned int i,j, nodesexamined,success;
unsigned int *isactive=cyclestore->isactive;
unsigned int *incomingedges=cyclestore->incomingedges;
unsigned int **graph=cyclestore->graph;

for(i=0;i<numnodes;i++){isactive[i]=1;} /** all nodes initially active */

/** copy the current network definition into graph[][]*/
for(i=0;i<numnodes;i++){for(j=0;j<numnodes;j++){graph[i][j]=dag->defn[i][j];}}

/** calc number of incoming edges for each child  and put into incomingedges**/
 get_numincomingedges(incomingedges,graph, numnodes);


/** find a node with no incoming edges - see lecture11.pdf in ../articles folder **/
nodesexamined=0;
success=1;
while(success){
        success=0;
	for(i=0;i<numnodes;i++){
        	  if(isactive[i] && !incomingedges[i]){/** if not yet examined and incoming edges */
                	isactive[i]=0; /** set this child to inactive*/
                	
                        /** remove all OUTGOING links from node i, e.g. wherever i is a parent */
                	for(j=0;j<numnodes;j++){graph[j][i]=0;}
                	get_numincomingedges(incomingedges,graph,numnodes);
           	        success=1; nodesexamined++;}
			}
           }
         

      if(nodesexamined==numnodes){return(0);/** no cycle */                               
      } else {/*Rprintf("=>%d %d\n",nodesexamined, numnodes);*/
              return(1);} /** found a cycle */
 

}

/** *************************************************************************************************/
/** v.small function but helps clarity in hascycle() */
/** *************************************************************************************************/
void get_numincomingedges(unsigned int *incomingedges,unsigned int **graph, unsigned int numnodes)
{ /** count up how many parents each child has **/
unsigned int i,j;
unsigned int numincomedge;
for(i=0;i<numnodes;i++){/** for each child */
        numincomedge=0;
 	for(j=0;j<numnodes;j++){numincomedge+=graph[i][j];
                                }
        incomingedges[i]=numincomedge;
        } 	

}

/** *************************************************************************************************/
/**  v.small function but helps clarity in main() */
/** *************************************************************************************************/
void calc_network_Score(storage *nodescore,network *dag,datamatrix *obsdata, double priordatapernode, int useK2, int verbose, SEXP R_labels)
{
 int i;
 int numnodes=dag->numNodes;
 double lognetworkscore=0.0;
 for(i=0;i<numnodes;i++){lognetworkscore+=calc_node_Score(nodescore,dag,obsdata,i,priordatapernode, useK2, verbose, R_labels);}
 dag->networkScore=lognetworkscore;
 
}

/** *************************************************************************************************/
/**  helps clarity in main() and checks to see whether the current node has previous been calculated */
/** *************************************************************************************************/
void calc_network_Score_reuse(storage *nodescore,network *dag,datamatrix *obsdata, double priordatapernode, int useK2, int verbose, SEXP R_labels,
                              struct database *prevNodes, int enforce_db_size)
{
 int i;
 int numnodes=dag->numNodes;
 double lognetworkscore=0.0;
 double curnodescore=0.0;
 for(i=0;i<numnodes;i++){
                         if(nodescoreisknown(dag,i,&curnodescore,prevNodes)){/** have previous calculated this node during search **/
                                                                 lognetworkscore+=curnodescore;
                         } else {/** not previously calculated so need to calculate this node */
                                 curnodescore=calc_node_Score(nodescore,dag,obsdata,i,priordatapernode, useK2, verbose, R_labels);
                                 /*Rprintf("nodescore=%f\n",curnodescore);*/ 
                                 storenodescore(dag,i,curnodescore,prevNodes, enforce_db_size); 
                                 lognetworkscore+=curnodescore;
                                 }
                         }
 dag->networkScore=lognetworkscore;
 
}


/** ***************************************************************************/
void init_network_score(storage *nodescore,network *dag){

int *parentindexes;
double *n_ij;
int *multipliers;
parentindexes=(int *)R_alloc(dag->maxparents,sizeof(int));
n_ij=(double *)R_alloc(dag->maxParentCombinations,sizeof(double));
multipliers=(int *)R_alloc(dag->maxparents,sizeof(int));

nodescore->parentindexes=parentindexes;
nodescore->n_ij=n_ij;
nodescore->multipliers=multipliers;/** this is from crossmultiply() - just setting up storage to avoid static allocation etc */

}
/** ***************************************************************************/
/** ***************************************************************************/
double calc_node_Score(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid,double priordatapernode, int useK2, int verbose, SEXP R_labels)
{
/** learn a single node - this code is esily validated against R lib deal**/

unsigned int totalparentcombs=0;
unsigned int i,j,k;
unsigned int numparents=0;
/*static int *parentindexes;*//** store the indexes in the network definition which are parents of the current node */
/*static int first=1;*/
int **parentcombs=dag->nodesparameters_lookup;/** will hold all parent combinations as indexes from 0,...k**/
double **parentcombsparams=dag->nodesparameters;/** parentcombinationindex->dirichlet param **/
int *levels=dag->numNodeLevels;/** will hold number of levels of each parent **/ 
int hit=0;
int response_obs=0;/** will hold the category of the observed value of the node at each data point**/
int response_obs_levels=obsdata->numVarlevels[nodeid];/** total number of categories in the response variable **/

double n_prior_ijk=0.0;
double n_prior_ij=0.0;
/*static double *n_ij;*//** this will hold the number of observations per parent combination **/ 
double logscore=0.0;
/*double score=0.0;*/

int *parentindexes=nodescore->parentindexes;/** just memory space **/
double *n_ij=nodescore->n_ij;/** just memory space **/


/** initialize **/
for(j=0;j<dag->maxparents;j++){levels[j]=0;} /** this will hold how many categories each individual parent node has **/


/** collect parents of this node and store their number of levels**/
for(j=0;j<dag->numNodes;j++){
              if(   dag->defn[nodeid][j]==1    /** got a parent so find how many levels each has **/
                 && numparents<dag->maxparents /** if numparents==dag->maxparents then we are done **/
                ){
								levels[numparents]=obsdata->numVarlevels[j];
								parentindexes[numparents++]=j;/** store index of parent **/
                                                                 
								}
				}

/** generate all possible combinations of levels of parents stored in parentcombs*/
totalparentcombs=crossproduct(nodescore,levels,parentcombs,numparents,dag->maxparents);

/*for(i=0;i<totalparentcombs;i++){for(j=0;j<numparents;j++){Rprintf("%d ",parentcombs[i][j]);}Rprintf("\n");}*/

/** as using the scratch space for each model need to zero the dirichlet parameters before any new learning **/
for(i=0;i<totalparentcombs;i++){	   n_ij[i]=0.0;/** reset to zero but only for nodes with parents */
	for(j=0;j<response_obs_levels;j++){parentcombsparams[i][j]=0.0;}
				} 

if(totalparentcombs==0){totalparentcombs=1;
                        parentcombs[0][0]=0;/** a dummy value to get the next bit to work not indep nodes**/
                        n_ij[0]=0.0;/** reset to zero for indep nodes **/
                        for(j=0;j<response_obs_levels;j++){parentcombsparams[0][j]=0.0;}/** reset parameters **/
			}

/** the pattern of levels we need to learn this node is in dag->nodesparameters_lookup[nodeorder][x][] **/
/** we loop over each array of dag->nodesparameters_lookup[nodeorder][x][] to update our parameters */

for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
	for(j=0;j<totalparentcombs;j++){/** try each combination of parents **/
		hit=0;/** used to check for a match **/
		for(k=0;k<numparents;k++){/** count up how many matches against observation and parent **/
				if(obsdata->data[i][parentindexes[k]]==parentcombs[j][k]+1){hit++;}		
				/** if value for this parent's value matches that in the current combination then increment counter**/
		}
		if(hit==numparents){/** have an exact match to the current data point and parent combination of levels**/
			/*printf("got an exact match to this point to dont need to check any more patterns\n");*/
			/** we now determine the observed value of the node i.e. the response category and save it **/
			response_obs=obsdata->data[i][nodeid];/** this is from 1,...maxnumcategories **/
                       
			/** the following two lines is the usual assumption that each data point has a weight of 1 **/
			parentcombsparams[j][response_obs-1]++;/** increment the appropriate dirichlet parameter **/
			n_ij[j]++;
			break;
					}
				}/** end of pattern check loop*/
		}/** end of data point loop*/


/** we now calculate the Bayesian Dirichlet Metric for the node **
    see p212 in machine learning 20,197-243 (1995) **/

/** we first use the uniformative parameter prior of priordatapernode split over all parent levels and response level **/
/** this is the BDeu metric of splitting the prior data equally across the node **/
if(useK2){
	n_prior_ijk=1.0;
        n_prior_ij=n_prior_ijk*response_obs_levels;
} else {
n_prior_ijk=priordatapernode/((double)(totalparentcombs*response_obs_levels));/** prior in dirichlet distn */
n_prior_ij=priordatapernode/(double)(totalparentcombs);/** prior per parent combination */
/*Rprintf("|%f %f|\n",n_prior_ijk,n_prior_ij);*/
}

/** print out all the combinations possible for the current node and the resulting frequencies **/
#ifdef DEBUG_LEARN_NODE
Rprintf("----------------------------------------------------------\n");
for(i=0;i<totalparentcombs;i++){/** for each parent combination **/
	Rprintf("\nNODE=%u\n",nodeid);
	Rprintf("total observations=%f,parentcombindex=%d numparents=%d\n",n_ij[i],i,numparents);
	for(k=0;k<numparents;k++){/** for each parent **/
		Rprintf("%u=%u ",parentindexes[k],parentcombs[i][k]+1);
	}
	Rprintf("| %u:",nodeid);
	for(k=0;k<response_obs_levels;k++){
	Rprintf("level %d = %f inc.prior=%f\t",k+1,parentcombsparams[i][k],(double)parentcombsparams[i][k]+n_prior_ijk);}
	Rprintf("\n");	
}
#endif
/** README **********************************************************************************************************/
/** a simple binomial response with beta prior has 
                                                  n_prior_ij  = alpha+beta
                                                  n_prior_ijk = alpha
                                                  n_prior_ijk = beta
                                                  n_ij        = n, total obs for that parentcomb
                                       parentcombsparams[j][k]= r
                                          response_obs_levels = 2,
hence doing it by hand we get a match to the 
                 (gamma(alpha+beta)/(gamma(alpha)gamma(beta))) x gamma(r+alpha)gamma(n-r+beta)/(gamma(n+alpha+beta))
**********************************************************************************************************************/

/** calulcate the logscore - use logs for numerical reasons as the gamma() will overflow **/


logscore=0.0;
for(j=0;j<totalparentcombs;j++){
	
	logscore+=lgammafn(n_prior_ij)-lgammafn(n_prior_ij+n_ij[j]);
	
	for(k=0;k<response_obs_levels;k++){
			logscore+=lgammafn(n_prior_ijk+parentcombsparams[j][k])-lgammafn(n_prior_ijk);			
	  
	}
	}
#ifdef DEBUG_LEARN_NODE
printf("n_prior_ij=%f n_prior_ijk=%f \nlogscore=%f %d\n",n_prior_ij,n_prior_ijk,logscore,totalparentcombs);
#endif	


if(verbose){/** we want to produce meaningful print outputs with correct labels and stuff **/

/*Rprintf("->%s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),9)));*/
Rprintf("----------------------------------------------------------\n");
Rprintf("------------PARAMETER ESTIMATES AND NODE SCORE------------\n");
Rprintf("----------------------------------------------------------\n");
for(i=0;i<totalparentcombs;i++){/** for each parent combination **/
	Rprintf("\n(CHILD) NODE = %s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),nodeid)));
	Rprintf("total observations=%.0f,parentcombindex=%d numparents=%d\n",n_ij[i],i,numparents);
        if(numparents>=1){Rprintf("PARENT NODES:\n");}
        for(k=0;k<numparents;k++){/** for each parent **/
		Rprintf("%s=%s ",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),parentindexes[k])),
                                 STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,parentindexes[k]+1),/** +1 needed see content of R_labels list() */
                                                                           parentcombs[i][k])));
	}
        Rprintf("\n");
	/*Rprintf("| %u:",nodeid);*/
        if(response_obs_levels<=2){Rprintf("Beta(level1,level2) parameters:\n");
        } else {Rprintf("Dirichlet(level1,level2,...) parameters:\n");}
	for(k=0;k<response_obs_levels;k++){
	Rprintf("level %d = %f inc.prior=%f\n",k+1,parentcombsparams[i][k],(double)parentcombsparams[i][k]+n_prior_ijk);}
	Rprintf("\n\n");	
}
        /** now for the node score **/
        Rprintf("n_prior_ij=%f n_prior_ijk=%f \nlogscore=%f\n",n_prior_ij,n_prior_ijk,logscore);

}




return(logscore);


}
/** ***************************************************************************/
/** ***************************************************************************/
int crossproduct(storage *nodescore,int *numlevels,int **dest, int n,int maxparents)
{
/** take vectors a[0,1], b[0,1,2] and produce a[0,1]xb[0,1,2] which is 
    00
    10
    01
    11
    02
    12
i.e. all possible combinations of level of a set of factors **/

/** numlevels is the number of categories in each vector being multiplied, dest is the
    place where all the combinations wil be stored, n is the number of vectors i.e. length of numlevels **/

/*static int *multipliers;
static int first=1;*/
int *multipliers=nodescore->multipliers;
int i,totalcombs,j,index;

/** only alloc memory the first time this fn is called **/
/*if(first){multipliers=(int *)malloc(maxparents*sizeof(int));first=0;}*/

multipliers[0]=1.0;
totalcombs=numlevels[0];
for(i=1;i<n;i++){multipliers[i]=multipliers[i-1]*numlevels[i-1];
                 totalcombs*=numlevels[i];}


for(i=0;i<totalcombs;i++){
                for(j=0;j<n;j++){
                          index=(i/multipliers[j])%(numlevels[j]);
			 /* printf("%d,",index);*/
			dest[i][j]=index;	}
			
			}



return(totalcombs);/** total number of combinations **/


}

/** *************************************************************************************/
/** read in a DAG definition from an R matrix and copy into C array in network struct **/
/** *************************************************************************************/
void setbanlist(network *dag,SEXP R_dag){

int i,j;

/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
/** dag->banlist memory allocated in build_init_dag() */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->banlist[i][j]=INTEGER(R_dag)[i+j*dag->numNodes];}} 

#ifdef DEBUG_dagbanlist   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->banlist[i][j]);} Rprintf("\n");}    
#endif

}
/** *************************************************************************************/
/** read in a DAG definition from an R matrix and copy into C array in network struct **/
/** *************************************************************************************/
void setretainlist(network *dag,SEXP R_dag){

int i,j;

/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
/** dag->banlist memory allocated in build_init_dag() */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->retainlist[i][j]=INTEGER(R_dag)[i+j*dag->numNodes];}} 

#ifdef DEBUG_dagretainlist   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->retainlist[i][j]);} Rprintf("\n");}    
#endif

}
/** *************************************************************************************/
/** read in a DAG definition from an R matrix and copy into C array in network struct **/
/** *************************************************************************************/
#ifdef OLD
void setstartlist(network *dag,SEXP R_dag){

int i,j;

/** this needs some work - copy a list of matrices INTEGER(R_dag)[i+j*dag->numNodes] becomes 
    INTEGER(VECTOR_ELT(R_dag,k)) for the kth matrix in the list**/

/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
/** dag->banlist memory allocated in build_init_dag() */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->startlist[i][j]=INTEGER(R_dag)[i+j*dag->numNodes];}} 

#ifdef DEBUG_dagstartlist   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->startlist[i][j]);} Rprintf("\n");}    
#endif

}
#endif
/** *************************************************************************************************/
/** *************************************************************************************************/
void init_random_dag(storage *nodescore,network *dag){

unsigned int maxlinkspossible=((dag->numNodes*dag->numNodes)-dag->numNodes);/** total number of links possible 
                                                                          removing child being also parent but ignoring cycles **/
unsigned int *ordertmp,i;
unsigned int **order=(unsigned int **)R_alloc( maxlinkspossible,sizeof(unsigned int*));
for(i=0;i<maxlinkspossible;i++){ordertmp=(unsigned int *)R_alloc( 2,sizeof(unsigned int)); order[i]=ordertmp;}
unsigned int *indexes=(unsigned int *)R_alloc( maxlinkspossible,sizeof(unsigned int));

nodescore->order=order;
nodescore->indexes=indexes;

}

/** *************************************************************************************************/
/** generate a random network from which to start a hill climbing search */
/** do this but always starting with a null network and randomly add arcs, checking for cycle each time **/
/** *************************************************************************************************/
void generate_random_dag(cycle *cyclestore,storage *nodescore,network *dag,unsigned nopermuts, unsigned int maxparents, SEXP R_shuffle, unsigned int offset,
                         unsigned int searchnum, SEXP R_dag_start)
{
 /** nopermuts is the number of links to add to the null network **/
/*static unsigned int **order,*ordertmp;
static unsigned int *indexes;*/
unsigned int **order=nodescore->order;
unsigned int *indexes=nodescore->indexes;
unsigned int i,j,row, numtried;
unsigned int nolinks=0;
unsigned int maxlinkspossible=((dag->numNodes*dag->numNodes)-dag->numNodes);/** total number of links possible 
                                                                          removing child being also parent but ignoring cycles **/
/*static int first=1;

if(first){order=(unsigned int **)malloc( maxlinkspossible*sizeof(unsigned int*));
	for(i=0;i<maxlinkspossible;i++){ordertmp=(unsigned int *)malloc( 2*sizeof(unsigned int)); order[i]=ordertmp;}
        indexes=(unsigned int *)malloc( maxlinkspossible*sizeof(unsigned int));*//** so can iterate over order[] **/
        /*first=0;}*/



nullnetworkdefn(dag,cyclestore, searchnum, R_dag_start);/** start off with a graph comprising of the given start graph */

row=0;
for(i=0;i<dag->numNodes;i++){
	for(j=0;j<dag->numNodes;j++){
		if(i!=j){order[row][0]=i;
                         order[row++][1]=j;
                         }
                                }
                         }

for(i=0;i<maxlinkspossible;i++){indexes[i]=INTEGER(R_shuffle)[i+offset];/*Rprintf("|%u-%u|",indexes[i],i);*/} 
/** roll shuffled R vector into this - offset is needed as all random indexes for all searches are in one vector  **/

numtried=0;/** this iterates through indexes[0]---[(maxlinkspossible-1)], indexes has been shuffled to make random arc choice */
while(nolinks<nopermuts && numtried<maxlinkspossible){/** keep trying until either have all links we want or else none to try */
        if(add_arc(cyclestore,dag,order[indexes[numtried]][0],order[indexes[numtried]][1],maxparents)){nolinks++;}
        numtried++;
}

if(nolinks<nopermuts){Rprintf("warning: not all init.permuts arcs can be added into initial network: arcs added=%u\n",nolinks);}
	
/*Rprintf("number of init.permuts arc=%d\n",nolinks);*/

}

/** **************************************************************************************************/
/** single step in hill-climbing search **************************************************************/
/** **************************************************************************************************/
void hill_climb_iter(storage *nodescore,cycle *cyclestore,network *dag_orig, network *dag_scratch,network *dag_opt1,
		     network *dag_opt2, network *dag_opt3, unsigned int maxparents, datamatrix *obsdata, double priordatapernode,
		     int useK2, int verbose, SEXP R_labels, struct database *prevNodes, int enforce_db_size)
{
 /** take a passed network and find network score for all the networks which add one, remove one, reverse one link **/
 /** iterate through each child and add a single link **/
/** reutrn the network with the best fit - don't compare this here with the current accepted network (dag_orig) */
unsigned int i,j;
double bestscore=-HUGE_VAL;
double curscore=-HUGE_VAL;
dag_opt1->networkScore=-HUGE_VAL;/** worst possible goodness of fit **/
dag_opt2->networkScore=-HUGE_VAL;
dag_opt3->networkScore=-HUGE_VAL;


/** 1. add an arc **/
 for(i=0;i<dag_orig->numNodes;i++){
	for(j=0;j<dag_orig->numNodes;j++){
			if(i!=j){/**child cannot also be its parent! **/ /*printf("\n%d %d ",i,j);*/
                           copynetworkdefn(dag_orig,dag_scratch); /** this could be streamlined since to re-copies full matrix **/        
                           if(add_arc(cyclestore,dag_scratch,i,j,maxparents)){ /** add an arc IF IT DOES NOT CREATE A CYCLE **/                        
                              /** adding arc was successful so get networkscore **/
                              calc_network_Score_reuse(nodescore,dag_scratch,obsdata,priordatapernode, useK2, 0, R_labels,prevNodes,enforce_db_size);
			      curscore=dag_scratch->networkScore;   
                              if(curscore>bestscore){bestscore=curscore;
                                                     copynetworkdefn(dag_scratch,dag_opt1);/** replace current graph with new best graph */
                                                     dag_opt1->networkScore=bestscore;
						     /*Rprintf("adding: got new best score %f\n", dag_opt1->networkScore);*/}
                              /*write_network(dag_scratch,outfilePtr1);*/
                              /* fprintf(outfilePtr1,"added arc %f\n\n",curscore);*/   
                             } /** end of i!=j if **/
                                 }
			}
		}

/** note - dont need to re-set bestscore since if add arc model is better then no point in doing anything anyway*/
/*bestscore=dag_orig->networkScore;*/
/** 2. remove an arc **/
 for(i=0;i<dag_orig->numNodes;i++){
	for(j=0;j<dag_orig->numNodes;j++){
			if(i!=j){/**child cannot also be its parent! **/
                           copynetworkdefn(dag_orig,dag_scratch);        
                           if(remove_arc(dag_scratch,i,j)){ /** remove arc **/                        
                              /** remove arc was successful so get networkscore **/
                             calc_network_Score_reuse(nodescore,dag_scratch,obsdata,priordatapernode, useK2, 0, R_labels,prevNodes,enforce_db_size);
			      curscore=dag_scratch->networkScore;
				if(curscore>bestscore){bestscore=curscore;
                                                     copynetworkdefn(dag_scratch,dag_opt2);/** replace current graph with new best graph */
							dag_opt2->networkScore=bestscore;
						     /*Rprintf("removing: got new best score %f\n", dag_opt2->networkScore);*/}
                              /*write_network(dag_scratch,outfilePtr1);*/
                               /*fprintf(outfilePtr1,"removed arc %f\n\n",curscore);*/   
                             } /** end of i!=j if **/
                                 }
			}
		}

/*bestscore=dag_orig->networkScore;*/
/** 3. reverse an arc **/
 for(i=0;i<dag_orig->numNodes;i++){
	for(j=0;j<dag_orig->numNodes;j++){
			if(i!=j){/**child cannot also be its parent! **/
                           copynetworkdefn(dag_orig,dag_scratch);         
                           if(reverse_arc(cyclestore,dag_scratch,i,j,maxparents)){ /** reverse an arc **/                        
                              /** reversing arc was successful so get networkscore **/
                              calc_network_Score_reuse(nodescore,dag_scratch,obsdata,priordatapernode, useK2, 0, R_labels,prevNodes,enforce_db_size);
			      curscore=dag_scratch->networkScore;
				if(curscore>bestscore){bestscore=curscore;
                                                     copynetworkdefn(dag_scratch,dag_opt3);/** replace current graph with new best graph */
							dag_opt3->networkScore=bestscore;
						     /*Rprintf("reversing: got new best score %f\n", dag_opt3->networkScore);*/} 
                              /*write_network(dag_scratch,outfilePtr1);*/
                               /*fprintf(outfilePtr1,"reversing arc %f\n\n",curscore);*/   
                             } /** end of i!=j if **/
                                 }
			}
		}


/** finally find out which of add (dag_opt1), remove (dag_opt2), reverse (dag_opt3) gave the best one-step model and then update dag_orig */
/** adding came up with the best single step model */
if(   dag_opt1->networkScore > dag_opt2->networkScore &&
      dag_opt1->networkScore > dag_opt3->networkScore){
                                                         copynetworkdefn(dag_opt1,dag_orig);
                                                         dag_orig->networkScore=dag_opt1->networkScore;
                                                         if(verbose){Rprintf("added arc...\n");}return;}
/** removing came up with the best single step model */
if(   dag_opt2->networkScore > dag_opt1->networkScore &&
      dag_opt2->networkScore > dag_opt3->networkScore){
                                                         copynetworkdefn(dag_opt2,dag_orig);
                                                         dag_orig->networkScore=dag_opt2->networkScore;
                                                         if(verbose){Rprintf("removed an arc...\n");}return;}
/** reversing came up with the best single step model */
if(   dag_opt3->networkScore > dag_opt1->networkScore &&
      dag_opt3->networkScore > dag_opt2->networkScore){
                                                         copynetworkdefn(dag_opt3,dag_orig);
                                                         dag_orig->networkScore=dag_opt3->networkScore;
                                                         if(verbose){Rprintf("reversed an arc...\n");}return;}
  /** at the end of this dag has the new best network in it **/
} 
/** *************************************************************************************************/
/** set initial DAG to the retain list - which must be provided even if its just a matrix of zeros  */
/** *************************************************************************************************/
void nullnetworkdefn(network *dag,cycle *cyclestore, unsigned int searchnum, SEXP R_dag_start)
{
 /** copy a matrix - dag - passed from R and set this each to dag->defn **/
 int i,j;
 /*for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){dag->defn[i][j]=dag->startlist[i][j];}}*/
 for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){dag->defn[i][j]=INTEGER(VECTOR_ELT(R_dag_start,searchnum))[i+j*dag->numNodes];}}   
 /*INTEGER(VECTOR_ELT(R_dag,k)) kth entry in a list*/
 /** check if this is acyclic **/
 if(hascycle(cyclestore,dag)){error("The model definition in matrix start.m contains a cycle - all graphs must be acyclic");}
}

/** *************************************************************************************************/
/** add an arc to an existing network */
/** *************************************************************************************************/
unsigned int add_arc(cycle *cyclestore,network *dag, unsigned int child, unsigned int parent, unsigned int maxparents)
{
 unsigned int numnodes=dag->numNodes;
 unsigned int i;
 unsigned int numparents=0;
 /** step 1. check whether arc is already present **/
 if(!dag->defn[child][parent] && !dag->banlist[child][parent]){/** this link is not yet present and NOT banned so add it **/
              for(i=0;i<numnodes;i++){numparents+=dag->defn[child][i];} /** check that child does not already have maxparents **/
              if(numparents<maxparents){dag->defn[child][parent]=1;} /** add link */	
              if(hascycle(cyclestore,dag)){/** cant add this arc - revert back to original graph **/
                                 /*printf(" hascycle\n");*/
                                      dag->defn[child][parent]=0;}
return(dag->defn[child][parent]); /** 1 if could add link, 0 if otherwise (due to either maxparents or cycle */
} else {return(0);}





}
/** *************************************************************************************************/
/** remove an arc from an existing network */
/** *************************************************************************************************/
 unsigned int remove_arc(network *dag, unsigned int child, unsigned int parent)
{
 if(dag->defn[child][parent] && !dag->retainlist[child][parent]){/** if arc currently in the graph and NOT in retain list**/
            dag->defn[child][parent]=0; /** set to zero **/
            return(1);/** successfully removed the arc **/
 } else {return(0);} /** no arc to remove **/

}

/** *************************************************************************************************/
/** reverse an arc in existing network */
/** *************************************************************************************************/
unsigned int reverse_arc(cycle *cyclestore,network *dag, unsigned int child, unsigned int parent, unsigned int maxparents)
{
 unsigned int numnodes=dag->numNodes;
 unsigned int i;
 unsigned int numparents=0;
 /** step 1. check whether arc is already present **/
 if(    dag->defn[child][parent] /** link is present */
     && !dag->defn[parent][child] /** reverse link is not present */
     && !dag->banlist[parent][child] /** reverse link is not banned */
     && !dag->retainlist[child][parent] /** current link is not in retain list */
    ){
              for(i=0;i<numnodes;i++){numparents+=dag->defn[parent][i];} /** check that reversed link will not already have maxparents **/
              if(numparents<maxparents){dag->defn[parent][child]=1;
                                        dag->defn[child][parent]=0;} /** reverse link */	
              if(hascycle(cyclestore,dag)){/** cant add this arc - revert back to original graph **/
                                      dag->defn[parent][child]=0;
                                      dag->defn[child][parent]=1;}
		
   return(dag->defn[parent][child]); /** 1 if could reverse link, 0 if otherwise */
  } else {return(0);}


}
/** *************************************************************************************************/
/** try and make this inline since its a v.small function but helps clarity COPY DAG to DAG_SCRATCH */
/** *************************************************************************************************/
void copynetworkdefn(network *dag,network *dag_scratch)
{
 int i,j;
 for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){dag_scratch->defn[i][j]=dag->defn[i][j];}}
 
}

