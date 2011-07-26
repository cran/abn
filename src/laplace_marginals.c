#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"        
#include "network.h"
#include "laplace.h"
#include "laplace_marginals.h"
/** gsl functions */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define PRINTGSL1
/** *************************************************************************************************/
/**  calculate marginal posterior distribution for a SINGLE variable in a given node                */
/**  returns results as a matrix first col x, second col pdf(x)                                     */
/** *************************************************************************************************/
void calc_network_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd,
				 int nodenum, int varnum, gsl_matrix *posterior)
{
 int i;
 int numnodes=dag->numNodes;
 double lognetworkscore=0.0;
 double tmp=0;
 /** nodenum is the index of the node - the individual glm model **/
 /** varnum is the index of the variable in the individual glm model whose density is needed **/
 /** get node score since this is needed as the denominator for the posterior density **/ 
 lognetworkscore=calc_node_Score_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix, R_labels, priormean, priorsd);
 /** now get the numerator - a vector **/
 calc_node_Marginals_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix, R_labels, priormean, priorsd, varnum, posterior, lognetworkscore);
                        
           
 }           
     
/** ****************************************************************************************************
 ***** calc an individual logistic regression model - see laplace.c for full margLike calc - also easier code but similar
 *******************************************************************************************************/
void calc_node_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd, int varnum, gsl_matrix *posterior,
				double denominator)
{
 int i,j,k,ss,status;
 int iter=0;
 int numparents=0;
 int *levels=dag->numNodeLevels;/** will hold number of levels of each parent **/ 
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 
 int *parentindexes=nodescore->parentindexes;/** just memory space **/
 
 gsl_vector *Y,*myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*dgvalue,*term1,*term2,*term3,*vectmp3long,*vecpriormean,*vecpriorsd,*betafull,*dgvaluesfull;
 gsl_matrix *hessgvalue,*mattmp1,*mattmp2,*datamatrix,*hessgvaluefull;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm;
 
 
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

  datamatrix=gsl_matrix_alloc(obsdata->numDataPts,numparents+1);
  designmatrix->datamatrix=datamatrix;
  Y=gsl_vector_alloc(obsdata->numDataPts);
  designmatrix->Y=Y;
	vecpriormean=gsl_vector_alloc(numparents+1);
  designmatrix->priormean=vecpriormean;	
	vecpriorsd=gsl_vector_alloc(numparents+1);
  designmatrix->priorsd=vecpriorsd;
 /** create design matrix - copy relevant cols from the observed data **/
 /** int** designmatrix is just used as storage space, fill up from left cols across until as far as needed */
 for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
   /*designmatrix->data[i][0]=1;*//** intercept term **/
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0);
   /*designmatrix->Y[i]=obsdata->data[i][nodeid]-1;*//** copy values at node - response values - into vector: NOTE: -1 is as raw is coded 1,2, not 0,1 */
   gsl_vector_set(designmatrix->Y,i,obsdata->data[i][nodeid]-1);
   for(k=0;k<numparents;k++){
                              /*designmatrix->data[i][k+1]=obsdata->data[i][parentindexes[k]]-1;*/
			      gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->data[i][parentindexes[k]]-1);
			      /** copy the relevant cols of the observed data into the design matrix -off set by first col for intercept
			          and the value of the cell is reduced by one since the raw is coded 1,2 not 0,1 **/
                            }
                         }   
                        
   designmatrix->numparams=numparents+1;/** +1 for intercept**/
   /** now set the priormean and priorsd vector - choose the correct prior values */
   /*designmatrix->priormean[0]=priormean[0];*/
   gsl_vector_set(designmatrix->priormean,0,priormean[0]);
   /*designmatrix->priorsd[0]=priorsd[0];*/
   gsl_vector_set(designmatrix->priorsd,0,priorsd[0]);

   for(k=0;k<designmatrix->numparams-1;k++){gsl_vector_set(designmatrix->priormean,k+1,priormean[parentindexes[k]+1]);/** +1 since first entry is constant */
                                            gsl_vector_set(designmatrix->priorsd,k+1,priorsd[parentindexes[k]+1]);/** +1 since first entry is constant */
  }

  /** down to here is as for the network score calc which is an integral across all parameters - we now adjust this to that its across all parameters
      minus one, where this one is fixed at values across a grid **/
  /** SPECIAL CASE if only a model with a single parameter then no integration required just evaluation of (-1/n)*g() **/
  /** ********************************************************************************************************************/
  switch(designmatrix->numparams){
    case 1:{/** only a constant term **/  
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    dgvalue = gsl_vector_alloc (designmatrix->numparams);/** will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
       
   /** now send */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.vectmp3long=vectmp3long;
   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;
   gparams.mattmp1=mattmp1;
   gparams.mattmp2=mattmp2;
   
   myBeta = gsl_vector_alloc (designmatrix->numparams);   
   n=obsdata->numDataPts;
   
   for(i=0;i<posterior->size1;i++){
      gsl_vector_set(myBeta,0,gsl_matrix_get(posterior,i,0));
      laplace_g(myBeta,&gparams, &gvalue);
      logscore= -n*gvalue;
      gsl_matrix_set(posterior,i,1,exp(logscore-denominator));
      R_CheckUserInterrupt();/** allow an interupt from R console */ 
      }
      
      gsl_vector_free(Y);
      gsl_vector_free(myBeta);
      gsl_vector_free(vectmp1);
      gsl_vector_free(vectmp2);
      gsl_vector_free(vectmp1long);
      gsl_vector_free(vectmp2long);
      gsl_vector_free(dgvalue);
      gsl_vector_free(term1);
      gsl_vector_free(term2);
      gsl_vector_free(term3);
      gsl_vector_free(vectmp3long);
      gsl_vector_free(vecpriormean);
      gsl_vector_free(vecpriorsd);
      gsl_matrix_free(mattmp1);
      gsl_matrix_free(mattmp2);
      gsl_matrix_free(datamatrix);
      
    break;  
    }
    /** ********************************************************************************************************************/
    /** THIS IS THE MAIN CASE AND REST OF FUNCTION IS IN HERE***************************************************************/
    default:{/** this is for all models which contain more than just a constant **/
  
  /** allocate only once here since same dimension within one node for the marginals**/
  /** GENERAL IDEA - keep the same dimensions as in the full margLik calc but drop off terms at the end if needed */
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    dgvalue = gsl_vector_alloc (designmatrix->numparams-1);/** IMPORTANT: -1 since now a marginal calculation will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    betafull = gsl_vector_alloc (designmatrix->numparams);/** this will hold the re-build full vector of all parameters */
    hessgvaluefull = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);/**  will hold hessian matrix **/
    dgvaluesfull = gsl_vector_alloc (designmatrix->numparams);/** IMPORTANT: -1 since now a marginal calculation will hold partial derivates **/
    
    n=obsdata->numDataPts;
    m=designmatrix->numparams-1;/** IMPORTANT: -1 since now a marginal calculation **/
    FDF.f = &laplace_dg_marg;
    FDF.df = &laplace_hessg_marg;
    FDF.fdf = &wrapper_fdf_marg;
    FDF.n = designmatrix->numparams-1;
    FDF.params = &gparams;
    myBeta = gsl_vector_alloc (designmatrix->numparams-1);/** this will hold the parameter point estimates */   
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams-1,designmatrix->numparams-1);/**  IMPORTANT: -1 since now a marginal calculation will hold hessian matrix **/
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams-1);
    perm = gsl_permutation_alloc (m);
    
   /** now send */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.vectmp3long=vectmp3long;
   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;
   gparams.mattmp1=mattmp1;
   gparams.mattmp2=mattmp2;
   gparams.betafull=betafull;
   gparams.dgvalues=dgvaluesfull;
   gparams.hessgvalues=hessgvaluefull;
   gparams.betafixed=0.0;/** these will be changed in loop below*/
   gparams.betaindex=varnum;/** this is fixed - the variable for which the posterior is calculated **/
 
   generate_inits(myBeta,designmatrix); /** generate initial estimates for the remaining variable - not the posterior variable **/  
   
   /** POSTERIOR DENSITY CALC STARTS HERE **/ 
   for(i=0;i<posterior->size1;i++){
    gparams.betafixed=gsl_matrix_get(posterior,i,0);/** this is passed from R - the fixed value for the posterior variable **/  
     gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
     iter=0; 
       do
         {
           iter++;
           status = gsl_multiroot_fdfsolver_iterate (s);
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, 1e-7);
         }
       while (status == GSL_CONTINUE && iter < 100);
     
       if(status != GSL_SUCCESS){Rprintf ("Zero finding error: status = %s\n", gsl_strerror (status));exit(1);}
      gsl_vector_memcpy(myBeta,s->x);
    /** we now have all the individual parts so put it together to the laplace approx */
    laplace_g_marg(myBeta,&gparams, &gvalue);
    laplace_hessg_marg(myBeta,&gparams, hessgvalue);
    gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
    logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2)*log((2*M_PI)/n); /** this is the final value */
    gsl_matrix_set(posterior,i,1,exp(logscore-denominator));
    /*Rprintf("%f %5.10f\n",gsl_matrix_get(posterior,i,0),gsl_matrix_get(posterior,i,1));*/
   R_CheckUserInterrupt();/** allow an interupt from R console */ 
   }
  
    /*** Last Step before return - free all the gsl vectors, matrices, other etc **/
   gsl_vector_free(Y);
   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(dgvalue);
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_vector_free(vecpriormean);
   gsl_vector_free(vecpriorsd);
   gsl_vector_free(betafull);
   gsl_vector_free(dgvaluesfull);
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp1);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(datamatrix);
   gsl_matrix_free(hessgvaluefull);
   gsl_permutation_free(perm);
   gsl_multiroot_fdfsolver_free (s);
     
    }} /** end of switch **/
   
}

/** ***************************************************************************************************************
******************************************************************************************************************* 
** laplace method = int^b_a exp(-lambda g(y)) h(y) dy = exp(-lambda g(y*)) h(y*) (2PI/lambda)^(d/2) det(hess)^(1/2)
** lambda = sample size n, g(y) = -(1/n)* log( f(D|betas)f(betas) ) e.g. -(1/n)* log (like*prior)
*******************************************************************************************************************
******************************************************************************************************************/

/** **************************************************************************************************************/
/** g(y) = -(1/n)* log( f(D|betas)f(betas) */ 
/** **************************************************************************************************************/
int laplace_g_marg (const gsl_vector *betashort, void *params,double *gvalue)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
       gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       double n=Y->size;/** no. observations **/
       double m=X->size2;
        /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       double term1=0;
       double term2=0;
       double term3=0;
       double storedbl1;
       int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
     
      /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }	 
     
     /** DO IN THREE PARTS - term1, term2, term3 */
     /** R code "term2<-sum( log(1/(sqrt(2*pi)*sd.loc)) );" **/
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
     
     /** Calc this in parts: R code "term3<- sum( (-1/(2*sd.loc*sd.loc))*(mybeta-mean.loc)*(mybeta-mean.loc) );" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,vectmp1);/** copy vectmp1 to vectmp2 **/
     gsl_vector_mul(vectmp2,vectmp1);/** square all elements in vectmp1 and store in vectmp2 */
     gsl_vector_memcpy(vectmp1,priorsd);
     gsl_vector_mul(vectmp1,priorsd);/** square all elements in priorsd and store in vectmp1 */
     gsl_vector_div(vectmp2,vectmp1);/** vectmp2/vectmp1 and store in vectmp2 **/
     gsl_vector_scale(vectmp2,-0.5); /** scale by -1/2 */
     gsl_vector_set_all(vectmp1,1.0); /** ones vector */
     gsl_blas_ddot (vectmp2, vectmp1, &term3);/** DOT product simply to calcu sum value */
     /*term3 = -(1.0/(2.0*priorsd*priorsd))*storedbl1;*/
          
     
     /** Rcode  Y%*%(X%*%mybeta)-sum(log(1+exp(X%*%mybeta)));  */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     term1+=storedbl1;
     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-log(1.0+exp(gsl_vector_get(vectmp1long,i))));
                                      } /** vectmp2 holds -log(1+exp(X%*%mybeta)) */
     gsl_vector_set_all(vectmp1long,1.0); /** ones vector */  
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** DOT product simply to calc -sum(log(1+exp(X%*%mybeta))) */
     term1+=storedbl1;
     
     *gvalue=(-1.0/n)*(term1+term2+term3);

       return GSL_SUCCESS;
     }
      
   
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_dg_marg (const gsl_vector *betashort, void *params, gsl_vector *dgvaluesshort)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_vector *term3 = ((struct fnparams *) params)->term3;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       double n=Y->size;/** no. observations **/

       int i=0; double tmp=0;int col;
       /** beta are the parameters values at which the function is to be evaluated **/
     /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}
       
       
     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }	 
         
     /** DO IN THREE PARTS - term1, term2, term3 */
      /** term3 (beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0); 
     gsl_vector_memcpy(term1,vectmp1);
  
     /** Rcode  -sum(log(1+exp(X%*%mybeta)));  */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-exp(gsl_vector_get(vectmp1long,i))/(1+exp(gsl_vector_get(vectmp1long,i))) );
                                      } /** vectmp2long holds exp(X%*%mybeta)/(1+exp(X%*%mybeta) */
     
     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp2long, 0.0, vectmp1);/** vectmp1long hold X%*%mybeta **/ 
     gsl_vector_memcpy(term2,vectmp1);
     
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);
     gsl_vector_memcpy(term3,vectmp1);
 
     /*Rprintf("==%f %f %f\n",gsl_vector_get(term1,0),gsl_vector_get(term2,0),gsl_vector_get(term3,0));*/
 
     gsl_vector_add(term1,term2);/** add term 2 to term 1 */
     gsl_vector_add(term1,term3);/** add term 3 to term 1 */
     gsl_vector_scale(term1,-1.0/n); 
     
     /** need to drop one cell in term1 before copying back */
     /** create an adjusted term1 which contains the term1 without the  re-inserted at the correct place **/
    col=0;
     for(i=0;i<beta->size;i++){
       if(i!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_vector_get(term1,i);
	 col=i;
	 if(i>betaindex){col=i-1;} 
                               gsl_vector_set(dgvaluesshort,col,tmp);}
	}
       
       
 return GSL_SUCCESS;
     }
       
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_hessg_marg (const gsl_vector *betashort, void *params, gsl_matrix *hessgvaluesshort)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
        gsl_vector *vectmp3long = ((struct fnparams *) params)->vectmp3long;/** numobs long **/
       /*const gsl_vector *priormean = ((struct fnparams *) params)->priormean;*/
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_matrix *hessgvalue = ((struct fnparams *) params)->hessgvalues;
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** no. params to estimate*/
       double tmp1=0;double tmp2=0;double tmp3=0;

       int i=0;int j=0;int k=0;int row,col;double tmp;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
     
       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}
       
       
     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }	 
     
     /** do in multiple parts - need to do element operations first */
     /** first exp(Xb) */ 
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     for(i=0;i<vectmp1long->size;i++){
           tmp1=exp(gsl_vector_get(vectmp1long,i));
           tmp2=1+exp(gsl_vector_get(vectmp1long,i));
           tmp3=exp(2.0*gsl_vector_get(vectmp1long,i));
            gsl_vector_set(vectmp2long,i,(tmp1*tmp2-tmp3)/(n*tmp2*tmp2));
                                      } /** vectmp2long holds the main complicated term */
     
     gsl_matrix_memcpy(mattmp1,X);/** make a copy of X*/
     gsl_matrix_mul_elements (mattmp1, X);/* calc X^2 is in mattmp1*/
     
     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, mattmp1, vectmp2long, 0.0, vectmp1);/** vecttmp2long hold Xij^2*complicated **/ 
     gsl_vector_memcpy(term1,vectmp1);
    
     
     gsl_vector_set_all(term2,0.0); /** zeros vector */
     
     gsl_vector_memcpy(vectmp1,priorsd);/** copy priorsd in vectmp1 **/
     gsl_vector_mul(vectmp1,priorsd);/** square priorsd */
     gsl_vector_scale(vectmp1,n);/** now have n*sigma^2 **/
     gsl_vector_set_all(vectmp2,1.0); /** ones vector */
     gsl_vector_div(vectmp2,vectmp1);/** get 1/(n*sigma^2) into vectmp2 **/
     gsl_vector_add(term2,vectmp2);/** add to term2*/
     
     gsl_vector_add(term1,term2);
     
     /*Rprintf("hess[1,1] at beta=%5.10f is %5.10f\n",gsl_vector_get(beta,0),gsl_vector_get(term1,0)); */
     
     
     /** STILL TO DO OFF DIAGONAL ELEMENTS - check for triangular?*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** dealt with j==k case above */
                          /** NOTE - vectmp2long is the SAME as in the j==k case so can use this directly **/
                          /** need X[,j]*X[,k] - element wise mult **/
                            gsl_matrix_get_col(vectmp1long,X,j); /** get col j in X **/
                            gsl_matrix_get_col(vectmp3long,X,k); /** get col k in X **/
                            gsl_vector_mul(vectmp1long,vectmp3long); /** element by element multiplication - result in vecttmp1long **/
                            
                            gsl_blas_ddot (vectmp1long, vectmp2long, gsl_matrix_ptr(hessgvalue,j,k));/** DOT product simply to calcu sum value */
                    } else {*gsl_matrix_ptr(hessgvalue,j,k)=gsl_vector_get(term1,j);}
                     }
                     }
                     
     /** need to drop a row and drop a col **/
     row=0;
     col=0;
     for(i=0;i<beta->size;i++){
        for(j=0;j<beta->size;j++){
       if(i!=betaindex && j!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_matrix_get(hessgvalue,i,j);
	 row=i;col=j;
	 if(i>betaindex){row=i-1;} 
	 if(j>betaindex){col=j-1;}
                               gsl_matrix_set(hessgvaluesshort,row,col,tmp);}
	}
       }
     
       
 return GSL_SUCCESS;
     }
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_fdf_marg (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_dg_marg(beta, gparams, dgvalues);
       laplace_hessg_marg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
