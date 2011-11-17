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
#include <gsl/gsl_sf_gamma.h>

#define PRINTGSL1
/** *************************************************************************************************/
/**  calculate marginal posterior distribution for a SINGLE variable in a given node                */
/**  returns results as a matrix first col x, second col pdf(x)                                     */
/** *************************************************************************************************/
void calc_network_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix,  const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale,
				 int nodenum, int varnum, int whichgaus,gsl_matrix *posterior,const int maxiters, const double epsabs)
{
 /*int numnodes=dag->numNodes;*/
 double lognetworkscore=0.0;
 /*double tmp=0;*/
 /** nodenum is the index of the node - the individual glm model **/
 /** varnum is the index of the variable in the individual glm model whose density is needed **/
 /** get node score since this is needed as the denominator for the posterior density **/ 
 int gaussiannodeid=whichgaus;
 /*Rprintf("whichgaus=%d\n",whichgaus);*/ 
 
 switch(obsdata->vartype[nodenum])  /** choose which type of node we have */
                         {
                           case 1:{ /** binary/categorical node */
                                    lognetworkscore=calc_node_Score_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,maxiters,epsabs);    
                                    calc_node_Marginals_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix,  priormean, priorsd,
                                                                priorgamshape,priorgamscale, varnum, posterior, lognetworkscore,maxiters,epsabs);                                                                    
                                     if(verbose){Rprintf("Binary node=%d score=%f\n", nodenum,lognetworkscore);}
                                     break;
                                   }
                         
                           case 0:{ /** gaussian node */
                                    lognetworkscore=calc_Gaussiannode_Score_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,gaussiannodeid,maxiters,epsabs);
                                    calc_Gaussiannode_Marginals_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale, varnum,  gaussiannodeid,posterior, lognetworkscore,maxiters,epsabs);                                                       
                                    if(verbose){Rprintf("Gaussian node=%d score=%f\n",nodenum,lognetworkscore);}
                                    break;
                                   }
                           default: {error("in default switch - should never get here!");}                                          
                         }  
 
 /*lognetworkscore=calc_node_Score_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix,  priormean, priorsd,priorgamshape,priorgamscale);*/
 

 /** now get the numerator - a vector **/

/* calc_node_Marginals_laplace(nodescore,dag,obsdata,nodenum,verbose, designmatrix,  priormean, priorsd,
                               priorgamshape,priorgamscale, varnum, posterior, lognetworkscore);
*/                        
           
 }           
     
/** ****************************************************************************************************
 ***** calc an individual logistic regression model - see laplace.c for full margLike calc - also easier code but similar
 *******************************************************************************************************/
void calc_node_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                         datamatrix *designmatrix, const double *priormean, const double *priorsd, const double *priorgamshape,
                          const double *priorgamscale, int varnum, gsl_matrix *posterior, double denominator,const int maxiters, const double epsabs)
{
 int i,j,k,ss,status;
 int iter=0;
 int numparents=0; 
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 
 int *parentindexes=nodescore->parentindexes;/** just memory space **/
 
 gsl_vector *Y,*myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*dgvalue,*term1,*term2,*term3,*vectmp3long,*vecpriormean,*vecpriorsd,*betafull/*,*dgvaluesfull*/;
 gsl_matrix *hessgvalue,*mattmp1,*mattmp2,*datamatrix,*hessgvaluefull,*mattmp3,*mattmp4;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm;
 gsl_permutation *initsperm;
 
 /** collect parents of this node and store their number of levels**/
 for(j=0;j<dag->numNodes;j++){
              if(   dag->defn[nodeid][j]==1    /** got a parent so find how many levels each has **/
                 && numparents<dag->maxparents /** if numparents==dag->maxparents then we are done **/
                ){
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
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0);
  /** copy values at node - response values - into vector: NOTE: -1 is as raw is coded 1,2, not 0,1 */
   gsl_vector_set(designmatrix->Y,i,obsdata->dataDouble[i][nodeid]-1);
   for(k=0;k<numparents;k++){/* now build design matric of explanatories */
	          
           switch(obsdata->vartype[parentindexes[k]])
	          {
	            case 1: {/** got discrete variable so map 1/2 to 0/1 using -1 operation*/
	                     gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]-1);break;}
	            case 0: {/** got gaussian - leave variable as is*/
	                     gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]); break;}
	            default: error("in default2 switch - should never get here!");          
            }
            
            /*gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]-1); */      
                            } /** end of explanatories **/
                         
    } /** end of data point loop */
                        
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
   /* dgvalue = gsl_vector_alloc (designmatrix->numparams-1);*//** IMPORTANT: -1 since now a marginal calculation will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    betafull = gsl_vector_alloc (designmatrix->numparams);/** this will hold the re-build full vector of all parameters */
    hessgvaluefull = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);/**  will hold hessian matrix **/
    /*dgvaluesfull = gsl_vector_alloc (designmatrix->numparams);*/ 
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
    
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
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.perm=initsperm;
   gparams.betafull=betafull;
   /*gparams.dgvalues=dgvaluesfull;*/
   gparams.hessgvalues=hessgvaluefull;
   gparams.betafixed=0.0;/** these will be changed in loop below*/
   gparams.betaindex=varnum;/** this is fixed - the variable for which the posterior is calculated **/
 
   /*generate_inits(myBeta,designmatrix);*/ /** generate initial estimates for the remaining variable - not the posterior variable **/  
   generate_inits_n(myBeta,&gparams);
    
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
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
     
       if(status != GSL_SUCCESS){Rprintf ("Zero finding error: status = %s at x=%f\n", gsl_strerror (status),gparams.betafixed);/*exit(1);*/}
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
  /* gsl_vector_free(dgvalue); */
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_vector_free(vecpriormean);
   gsl_vector_free(vecpriorsd);
   gsl_vector_free(betafull);
   /*gsl_vector_free(dgvaluesfull); */
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp1);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_matrix_free(datamatrix);
   gsl_matrix_free(hessgvaluefull);
   gsl_permutation_free(perm);
   gsl_permutation_free(initsperm);
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
/** GAUSSIAN NODE ALL BELOW HERE ***********************************************************************
********************************************************************************************************/   
/** ****************************************************************************************************
 *******************************************************************************************************/
void calc_Gaussiannode_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                         datamatrix *designmatrix, const double *priormean, const double *priorsd, const double *priorgamshape,
                          const double *priorgamscale, int varnum, int gaussiannodeid, gsl_matrix *posterior, double denominator,const int maxiters, const double epsabs)
{   
int i,j,k,ss,status;
 int iter=0;
 int numparents=0;
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 
 int *parentindexes=nodescore->parentindexes;/** just memory space **/
 
 gsl_vector *Y,*myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,/**dgvalue,*/*term1,*term2,*term3,*vectmp3long,*vecpriormean,*vecpriorsd,
            *vecpriorgamshape,*vecpriorgamscale,*localbeta,*betafull,*dgvaluesfull;
 gsl_matrix *hessgvalue,/**mattmp1,*/*mattmp2,*mattmp3,*mattmp4,*datamatrix,*hessgvaluefull;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm=0;
 
 /** collect parents of this node **/

 for(j=0;j<dag->numNodes;j++){
              if(   dag->defn[nodeid][j]==1    /** got a parent so find how many levels each has **/
                 && numparents<dag->maxparents /** if numparents==dag->maxparents then we are done **/
                ){
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
  vecpriorgamshape=gsl_vector_alloc(1); /** only 1 of these per node */
  designmatrix->priorgamshape=vecpriorgamshape;
  vecpriorgamscale=gsl_vector_alloc(1); /** only 1 of these per node */
  designmatrix->priorgamscale=vecpriorgamscale;
 
 /** create design matrix - copy relevant cols from the observed data **/
 /** int** designmatrix is just used as storage space, fill up from left cols across until as far as needed */
 for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
   /*designmatrix->data[i][0]=1;*//** intercept term **/
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0);
  /** copy values at node - response values - into vector */
  gsl_vector_set(designmatrix->Y,i,obsdata->dataDouble[i][nodeid]);
  
  
  for(k=0;k<numparents;k++){/* now build design matrice of explanatories */
	          
           switch(obsdata->vartype[parentindexes[k]])
	          {
	            case 1: {/** got discrete variable so map 1/2 to 0/1 using -1 operation*/
	                     gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]-1);break;}
	            case 0: {/** got gaussian - leave variable as is*/
	                     gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]); break;}
	            default: error("in default2b switch - should never get here!");          
            }
            
            /*gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->dataDouble[i][parentindexes[k]]-1); */      
                            } /** end of explanatories **/
                         
   } /** end of data point loop */   
                        
   designmatrix->numparams=numparents+1;/** +1 for intercept NOT including precision term **/
   /** now set the priormean and priorsd vector - choose the correct prior values */
   gsl_vector_set(designmatrix->priormean,0,priormean[0]);
   gsl_vector_set(designmatrix->priorsd,0,priorsd[0]);  
   for(k=0;k<designmatrix->numparams-1;k++){gsl_vector_set(designmatrix->priormean,k+1,priormean[parentindexes[k]+1]);/** +1 since first entry is constant */
                                            gsl_vector_set(designmatrix->priorsd,k+1,priorsd[parentindexes[k]+1]);/** +1 since first entry is constant */
                                            
   }
   
   gsl_vector_set(designmatrix->priorgamshape,0,priorgamshape[gaussiannodeid]); /** only one entry needed */
   gsl_vector_set(designmatrix->priorgamscale,0,priorgamscale[gaussiannodeid]); /** only one entry needed */
   /*Rprintf("index=%d shapeprior=%f scaleprior=%f\n",gaussiannodeid,gsl_vector_get(designmatrix->priorgamshape,0),
                                                                   gsl_vector_get(designmatrix->priorgamscale,0)); */  
   

   /** down to here is as for the network score calc which is an integral across all parameters - we now adjust this to that its across all parameters
      minus one, where this one is fixed at values across a grid **/
  /** ALWAYS HAVE AT LEAST TWO PARAMS: mean+precision if only a model with a single parameter then no integration required just evaluation of (-1/n)*g() **/
  /** ********************************************************************************************************************/
  switch(designmatrix->numparams+1){
    case 1:{ error("must always have at least two parameters - a mean term and a precision/variance!\n");break;}

   /** ********************************************************************************************************************/
    /** THIS IS THE MAIN CASE AND REST OF FUNCTION IS IN HERE***************************************************************/
    default:{/** **/
    /** allocate only once here since same dimension within one node for the marginals**/
   /** GENERAL IDEA - keep the same dimensions as in the full margLik calc but drop off terms at the end if needed */  
  
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    /*dgvalue = gsl_vector_alloc (designmatrix->numparams);*//** will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams+1-1,designmatrix->numparams+1-1);/** will hold hessian matrix NOTE: NOT used in solver
                                                                                            as it allocates its own storage based on the problem dimension**/
    /*mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,obsdata->numDataPts);*/
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions **/
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
    
    betafull = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the re-build full vector of all parameters */
    dgvaluesfull = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the re-build full vector of all parameters derivs */
    hessgvaluefull = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);/**  will hold hessian matrix **/

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
   gparams.priorgamshape=designmatrix->priorgamshape;
   gparams.priorgamscale =designmatrix->priorgamscale;
   /*gparams.mattmp1=mattmp1;*/
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.beta=localbeta;
   gparams.perm=initsperm;
   gparams.betafull=betafull;
   gparams.dgvaluesfull=dgvaluesfull;
   gparams.hessgvalues=hessgvaluefull;
   gparams.betafixed=0.0;/** these will be changed in loop below*/
   gparams.betaindex=varnum;/** this is fixed - the variable for which the posterior is calculated **/

    iter=0;       
    FDF.f = &laplace_gaus_dg_marg;
    FDF.df = &laplace_gaus_hessg_marg;
    FDF.fdf = &wrapper_gaus_fdf_marg;
    FDF.n = designmatrix->numparams+1-1;/** +1 for the precision term -1 as we drop off a single term as marginal post. estimation*/
    FDF.params = &gparams;  
    myBeta = gsl_vector_alloc (designmatrix->numparams+1-1);/** this will hold the parameter point estimates INC tau-precision -1=marginal*/
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams+1-1); /** +1 for the precision term -1=marginal */
    
    /*gparams.betafixed=gsl_matrix_get(posterior,0,0);
    generate_gaus_inits_marg(myBeta,&gparams); 
    laplace_gaus_g_marg(myBeta,&gparams, &gvalue);
    laplace_gaus_dg_marg (myBeta, &gparams, term1);
    laplace_gaus_hessg_marg(myBeta,&gparams, hessgvalue);
    */
     
   /** POSTERIOR DENSITY CALC STARTS HERE **/ 
   for(i=0;i<posterior->size1;i++){
    gparams.betafixed=gsl_matrix_get(posterior,i,0);/** this is passed from R - the fixed value for the posterior variable **/ 
    generate_gaus_inits_marg(myBeta,&gparams);/** generate initial estimates for the remaining variable - not the posterior variable **/
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
    
    iter=0; 
       do
         {
           iter++;
           status = gsl_multiroot_fdfsolver_iterate (s);
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
     
       if(status != GSL_SUCCESS){Rprintf ("Zero finding error: status = %s at x=%f\n", gsl_strerror (status),gparams.betafixed);/*exit(1);*/}
      gsl_vector_memcpy(myBeta,s->x);
    /** we now have all the individual parts so put it together to the laplace approx */
   
    laplace_gaus_g_marg(myBeta,&gparams, &gvalue);
    laplace_gaus_hessg_marg(myBeta,&gparams, hessgvalue);
    n=obsdata->numDataPts;
    m=designmatrix->numparams+1-1;/** -1 since a marginal calc. one less dimension in the integral **/
    perm = gsl_permutation_alloc (m);
    gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
    logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2)*log((2*M_PI)/n); /** this is the final value */
    gsl_matrix_set(posterior,i,1,exp(logscore-denominator));
    /*Rprintf("%f %5.10f\n",gsl_matrix_get(posterior,i,0),gsl_matrix_get(posterior,i,1));*/
   R_CheckUserInterrupt();/** allow an interupt from R console */ 
   
  }
   /*** Last Step before return - free all the gsl vectors, matrices, other etc **/
   /*Rprintf("IMPORTANT:still need to code up free to stop memory leaks\n");*/
   
   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(dgvaluesfull); 
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_vector_free(betafull);
   gsl_vector_free(localbeta);
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_matrix_free(hessgvaluefull);
   gsl_permutation_free(initsperm);
   gsl_permutation_free(perm);
   gsl_multiroot_fdfsolver_free (s);
   
   }
    } /** end of switch **/

   gsl_vector_free(Y);
   gsl_vector_free(vecpriormean);
   gsl_vector_free(vecpriorsd);
   gsl_vector_free(vecpriorgamshape);
   gsl_vector_free(vecpriorgamscale);
   gsl_matrix_free(datamatrix);

   
}   
/****************************************************************************************************************/ 
/**** Gaussian marginal functions down here *********************************************************************/
/****************************************************************************************************************/
/** **************************************************************************************************************/
int laplace_gaus_g_marg (const gsl_vector *betashort, void *params,double *gvalue)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
	    /* gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long; */
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;/** only long enough for beta terms - no precision */
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double term2=0;
       double term3=0;
       double term4=0;
       double term5=0;
       double term6=0;
       double storedbl1,storedbl2,storedbl3;
       double tau;      
       int i=0;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
            
        /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */
      
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/}
     
     /** same as in logistic model */
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
    /* Rprintf("term2 (Rterm3)=%f\n",term2);*/
     
     /** same as in logistic model */
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
    /* Rprintf("term3 (Rterm4)=%f\n",term3);*/     
     
     /** Need -0.5tau*(Y%*%Y+ (X%*%myBeta)%*%(X%*%myBeta)-2Y%*%X%*%myBeta) */
     
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     storedbl1= -2.0*storedbl1;/** now gives 2Y%*%X%*%myBeta */
     gsl_blas_ddot (vectmp1long, vectmp1long, &storedbl2);/** storebdl2 is (X%*%myBeta)%*%(X%*%myBeta) **/
     gsl_blas_ddot (Y, Y, &storedbl3); /** storebdl3 is Y%*%Y**/
     term4= -(tau/2.0)*(storedbl1+storedbl2+storedbl3);
     /*Rprintf("term4 (Rterm2)=%f\n",term4);*/
     
     term5= (n/2.0)*log(tau/(2.0*M_PI));
     /*Rprintf("term5 (Rterm1)=%f\n",term5);*/
     term6=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau)
	     -(tau/gsl_vector_get(priorgamscale,0));
    /* Rprintf("term6 (Rterm5)=%5.10f\n",term6);*/
     *gvalue=(-1.0/n)*(term2+term3+term4+term5+term6);

     /*Rprintf("loglike=%f\n",term5+term4);*/

     
     return GSL_SUCCESS;
     }
        
     
/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_dg_marg (const gsl_vector *betashort, void *params, gsl_vector *dgvaluesshort)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       /*gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
      /* gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1; */
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double tau;/*=gsl_vector_get(betaincTau,m);*//** extract the tau-precision from *beta */
       int i=0;int col;double tmp;
       double storedbl1/*,storedbl2,storedbl3,tmptau*/;
      /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_vector *dgvaluesfull = ((struct fnparams *) params)->dgvaluesfull;
       
       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */
      
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/} 
       
      /** do dg_beta terms first - store these in term1 vec*/ 
      /** -(beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0); 
     gsl_vector_memcpy(term1,vectmp1);/** the prior term in dg_dbeta **/
  
     /**  tau *sum(y_i X_ij - X_i beta X_ij) */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_scale(vectmp1,tau);/** mult by tau **/
     gsl_vector_memcpy(term2,vectmp1);/** the remaining term in dg_beta **/
     /** put the two parts of the dg_dbeta terms together */
     gsl_vector_add(term1,term2);
     gsl_vector_scale(term1,(-1.0/n));
     /** now store the dg_dbeta terms */
     for(i=0;i<m;i++){gsl_vector_set(dgvaluesfull,i,gsl_vector_get(term1,i));}
     
     /** dg_dtau **/
     gsl_vector_scale(vectmp1long,-1.0);/** from above vectmp1long is X*beta so this is -X*beta */
     gsl_vector_add(vectmp1long,Y);/** Y-X*beta **/
     gsl_vector_memcpy(vectmp2long,vectmp1long);
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** get sum((Y_i-X_i*beta)^2) using dot product */
     storedbl1= storedbl1*(-0.5);
     gsl_vector_set(dgvaluesfull,m,(-1.0/n)*(
                                          (n/(2.0*tau))+
				          storedbl1 + 
				          (gsl_vector_get(priorgamshape,0)-1.0)/tau
				          - (1.0/gsl_vector_get(priorgamscale,0)))
					  
                                          ); 
     
     /** need to drop one cell in term1 before copying back */
     /** create an adjusted term1 which contains the term1 without the  re-inserted at the correct place **/
    col=0;
     for(i=0;i<betafull->size;i++){
       if(i!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_vector_get(dgvaluesfull,i);
	 col=i;
	 if(i>betaindex){col=i-1;} 
                               gsl_vector_set(dgvaluesshort,col,tmp);}
	}
 
 /*Rprintf("end gaus_dg\n");
 for(i=0;i<betafull->size-1;i++){Rprintf("=%f\n",gsl_vector_get(dgvaluesshort,i));}	
 */
   return GSL_SUCCESS;
     }


/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_hessg_marg (const gsl_vector *betashort, void *params, gsl_matrix *hessgvaluesshort)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
      /* gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*//** numobs long **/
      /* const gsl_vector *priormean = ((struct fnparams *) params)->priormean; */
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
      /* const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale; */
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       /*gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
       gsl_matrix *mattmp3 = ((struct fnparams *) params)->mattmp3;
       gsl_matrix *mattmp2 = ((struct fnparams *) params)->mattmp2;
       int n=Y->size;/** no. observations **/
       int m=X->size2;/** number of coefficients excluding tau-precision */
       double tau;
       /*double tau=gsl_vector_get(betaincTau,m);*//** extract the tau-precision from *beta */
       int i=0;int row,col;double tmp;
       int j,k;
       double tmp1;
      /* double storedbl1,storedbl2,storedbl3,tmptau; */
      /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_matrix *hessgvaluefull = ((struct fnparams *) params)->hessgvalues;/** will hold "full hessian matrix" **/
       
       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */       
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/} 
       
  /** matrix X%*%t(X) - this is symmetrical so dg_bj_b_k=dg_b_k_b_j **/
  /*Rprintf("%d %d %d %d hesssize=%d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,hessgvalues->size1,hessgvalues->size2);*/
  gsl_matrix_memcpy(mattmp2,X);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                       1.0, X, mattmp2,
                       0.0, mattmp3);
  
  /*for(i=0;i<100;i++){Rprintf("%f \n",gsl_matrix_get(X,i,0));}*/
  
   /** matrix of hess excluding tau derivatives for the moment - note this is symmetrical*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** second derivatives */
                      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k));
		      *gsl_matrix_ptr(hessgvaluefull,j,k)=tmp1;
                    } else {
		      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k)-1.0/(gsl_vector_get(priorsd,j)*gsl_vector_get(priorsd,j)));
		      *gsl_matrix_ptr(hessgvaluefull,j,k)=tmp1;}
                     }
                     }
                    
  /** now for dg_dtau second deriv **/
  tmp1=(-1.0/n)*( (-n/(2.0*tau*tau)) - ( (gsl_vector_get(priorgamshape,0)-1)/(tau*tau)) );
  *gsl_matrix_ptr(hessgvaluefull,m,m)=tmp1;
  /*Rprintf("final cell is %d %d %f\n",m,m,*gsl_matrix_ptr(hessgvaluefull,m,m));*/
  
  /** now for dg_dtau_dbeta sum(y_i X_ij - X_i beta X_ij) - last ROW of hessian */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_memcpy(term1,vectmp1);/** the remaining term in dg_beta **/
     gsl_vector_scale(term1,(-1.0/n));
     
     /** last row in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvaluefull,m,j)=gsl_vector_get(term1,j);}
 
 /** dg_dbeta_dtau is same as dg_dtau_dbeta - note symmetry here **/
     
     /** last col in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvaluefull,j,m)=gsl_vector_get(term1,j);}
     
 /** now remove the row and col for the fixed variable */
 /** need to drop a row and drop a col **/
     row=0;
     col=0;
     for(i=0;i<betafull->size;i++){
        for(j=0;j<betafull->size;j++){
       if(i!=betaindex && j!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_matrix_get(hessgvaluefull,i,j);
	 row=i;col=j;
	 if(i>betaindex){row=i-1;} 
	 if(j>betaindex){col=j-1;}
                               gsl_matrix_set(hessgvaluesshort,row,col,tmp);}
	}
       }
 /*Rprintf("end gaus_dgHESS\n");
 for(i=0;i<betafull->size-1;i++){Rprintf("=%f\n",gsl_matrix_get(hessgvaluesshort,i,i));}*/	
 
  return GSL_SUCCESS;
} 
 
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int generate_gaus_inits_marg(gsl_vector *betashort,struct fnparams *gparams){
    
    /** betashort is all the beta's minus the one to be marginalised */
    /** general idea: find all the estimates and then simply drop the one that is fixed */    
    
    /** use least squares estimates to get started **/
    /** beta_hat= (X^T X)^{-1} X^T y **/
    
    const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_vector *vectmp1long = gparams->vectmp1long;/** numobs long **/
       gsl_vector *vectmp2long = gparams->vectmp2long;/** numobs long **/
       /*const gsl_vector *priormean = gparams->priormean;
       const gsl_vector *priorsd   = gparams->priorsd;
       const gsl_vector *priorgamshape = gparams->priorgamshape;
       const gsl_vector *priorgamscale = gparams->priorgamscale;*/
       /*gsl_vector *betalocal = gparams->beta;*//** this is beta without precision term*/
      /* gsl_vector *term1 = gparams->term1;
       gsl_vector *term2 = gparams->term2;
       gsl_vector *term3 = gparams->term3;*/
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_permutation *perm = gparams->perm;
      double n=Y->size;/** no. observations **/
     double m=X->size2;/** number of coefficients excluding tau-precision */
     int ss,i,j;
     double variance=0.0;
     gsl_vector *beta = gparams->betafull;/** will hold "full beta vector inc precision" **/   
     int betaindex = gparams->betaindex;
     
    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2);*/ 
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */                   
    gsl_linalg_LU_decomp(mattmp3,perm,&ss);
    gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */                   
    
    gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1); /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2); 
    
   for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(vectmp2,i));} /** set to Least squares estimate */
   
    /** now for variance estimate */
    /** first get y_hat estimate */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp2, 0.0, vectmp1long); /** vectmp1 is y_hat */ 
    gsl_vector_scale(vectmp1long,-1.0);/** - y_hat */
    gsl_vector_add(vectmp1long,Y);/** now have Y-y_hat (or -y_hat + Y) */
    gsl_vector_memcpy(vectmp2long,vectmp1long);
    gsl_blas_ddot (vectmp1long, vectmp2long, &variance);/** got sum((Y-Y_hat)^2) */
    variance=variance/(n-m);/** unbiased estimator using denominator n-#term in regression equation **/
    
    gsl_vector_set(beta,beta->size-1,1.0/variance);/** as are using precision parameterization not variance **/ 
    
    /** now copy into smaller vector and drop the variable which is fixed **/
    i=0;
    for(j=0;j<beta->size;j++){
      if(j!=betaindex){gsl_vector_set(betashort,i,gsl_vector_get(beta,j)); i++;}
    } /** set to Least squares estimate */
    
    /*Rprintf("passed varnum %d\n",betaindex);
    for(i=0;i<beta->size;i++){Rprintf("beta=%f\n",gsl_vector_get(beta,i));}
    for(i=0;i<betashort->size;i++){Rprintf("betashort=%f\n",gsl_vector_get(betashort,i));}
    */
    return GSL_SUCCESS;
}   
 

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_gaus_fdf_marg (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     { 
       laplace_gaus_dg_marg(beta, gparams, dgvalues);
       laplace_gaus_hessg_marg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }     
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
