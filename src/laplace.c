#include <R.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"        
#include "network.h"
#include "laplace.h"
#include "scorereuse.h"
/** gsl functions */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#define PRINTGSLa
#define IND

/** *************************************************************************************************/
/**  as above but score is calculated by evaluating the integrals directly rather than using conjugacy */
/** *************************************************************************************************/
void calc_network_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix, const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale,
         const int maxiters, const double epsabs, const int errverbose)
{
 int i;
 int numnodes=dag->numNodes;
 double lognetworkscore=0.0;
 double indnodescore=0.0;
 int gaussiannodeid=-1;/** this will hold the index of gaussian nodes e.g. 0,1,up to number of gaussian nodes -1 */
 /*Rprintf("Note: CONSTRAINED TO A SINGLE NODE\n");
 i=1;*/
 for(i=0;i<numnodes;i++){ 
                         switch(obsdata->vartype[i])  /** choose which type of node we have */
                         {
                           case 1:{ /** binary/categorical node */
                                    indnodescore=calc_node_Score_laplace(nodescore,dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,maxiters,epsabs); 
                                    if(verbose){Rprintf("Binary node=%d score=%f\n", i,indnodescore);}                                                                                         
                                    break;
                                   }
                         
                           case 0:{ /** gaussian node */
                                    indnodescore=0;
                                    gaussiannodeid++;
                                    indnodescore=calc_Gaussiannode_Score_laplace(nodescore,dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,gaussiannodeid,maxiters,epsabs);
                                    if(verbose){Rprintf("Gaussian node=%d score=%f\n",i,indnodescore);}
                                    break;
                                   }
                           default: {error("in default switch - should never get here!");}                                          
                         }                                     
                         R_CheckUserInterrupt();/** allow an interupt from R console */ 
                         /*if(verbose){Rprintf("individual node score=%f\n",indnodescore);}*/
                         lognetworkscore+=indnodescore;
                         } 
                         
 dag->networkScore=lognetworkscore;
 
 if(verbose){Rprintf("\n   #################################################################\n");
               Rprintf("   ###      log marginal likelihood for Model: %5.10f\n",lognetworkscore);
	             Rprintf("   #################################################################\n");
            }
            
 }           
/** *************************************************************************************************/
/**  as above but score is calculated by evaluating the integrals directly rather than using conjugacy */
/** *************************************************************************************************/
void calc_network_Score_laplace_reuse(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix, const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale,
         struct database *prevNodes,const int maxiters, const double epsabs, const int errverbose, int enforce_db_size)
{
 int i;
 int numnodes=dag->numNodes;
 double lognetworkscore=0.0;
 double indnodescore=0.0;
 int gaussiannodeid=-1;/** this will hold the index of gaussian nodes e.g. 0,1,up to number of gaussian nodes -1 */
 /*Rprintf("Note: CONSTRAINED TO A SINGLE NODE\n");
 i=1;*/
 for(i=0;i<numnodes;i++){ 
                         switch(obsdata->vartype[i])  /** choose which type of node we have */
                         {
                           case 1:{ /** binary/categorical node */
                                    if(nodescoreisknown(dag,i,&indnodescore,prevNodes)){/** have previous calculated this node during search **/
                                                                 /*lognetworkscore+=indnodescore;*/
                                    } else {/** not previously calculated so need to calculate this node */
                                    indnodescore=calc_node_Score_laplace(nodescore,dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,maxiters,epsabs);
                                    storenodescore(dag,i,indnodescore,prevNodes,enforce_db_size); }                                      
                                    if(verbose){Rprintf("Binary node=%d score=%f\n", i,indnodescore);}                                                                                         
                                    break;
                                   }
                         
                           case 0:{ /** gaussian node */
                                    gaussiannodeid++;
                                    indnodescore=0.0;
                                    if(nodescoreisknown(dag,i,&indnodescore,prevNodes)){/** have previous calculated this node during search **/
                                                                 /*lognetworkscore+=indnodescore;*/
                                    } else {/** not previously calculated so need to calculate this node */
                                    indnodescore=calc_Gaussiannode_Score_laplace(nodescore,dag,obsdata,i,errverbose, designmatrix, priormean, priorsd,
                                                                         priorgamshape,priorgamscale,gaussiannodeid,maxiters,epsabs);
                                    storenodescore(dag,i,indnodescore,prevNodes,enforce_db_size); }                                      
                                    if(verbose){Rprintf("Gaussian node=%d score=%f\n",i,indnodescore);}
                                    break;
                                   }
                           default: {error("in default switch - should never get here!");}                                          
                         }                                     
                         R_CheckUserInterrupt();/** allow an interupt from R console */ 
                         /*if(verbose){Rprintf("individual node score=%f\n",indnodescore);}*/
			 
                         lognetworkscore+=indnodescore;
                         }
                         
 dag->networkScore=lognetworkscore;
 
 if(verbose){Rprintf("\n   #################################################################\n");
               Rprintf("   ###      log marginal likelihood for Model: %5.10f\n",lognetworkscore);
	             Rprintf("   #################################################################\n");
            }
            
 }     
/** ****************************************************************************************************
 ***** calc an individual logistic regression model 
 *******************************************************************************************************/
double calc_node_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix, const double *priormean, const double *priorsd,
                                 const double *priorgamshape, const double *priorgamscale,const int maxiters, const double epsabs)
{
 int i,ss,status,status2,status_inits;
 unsigned int k,j;
 int iter=0;
 unsigned int numparents=0;
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 int *parentindexes=nodescore->parentindexes;/** just memory space **/
 
 gsl_vector *Y,*myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*dgvalue,*term1,*term2,*term3,*vectmp3long,*vecpriormean,*vecpriorsd;
 gsl_matrix *hessgvalue,*mattmp1,*mattmp2,*datamatrix,*mattmp3,*mattmp4;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 
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
 /** create design matrix - copy relevant cols from the observed data **/
 /** int** designmatrix is just used as storage space, fill up from left cols across until as far as needed */
 for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
   /*designmatrix->data[i][0]=1;*//** intercept term **/
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0);
  /** copy values at node - response values - into vector: NOTE: -1 needed as raw is coded 1,2, not 0,1 */
  gsl_vector_set(designmatrix->Y,i,obsdata->dataDouble[i][nodeid]-1);
   
   for(k=0;k<numparents;k++){/* now build design matrice of explanatories */
	          
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
   gsl_vector_set(designmatrix->priormean,0,priormean[0]);
   gsl_vector_set(designmatrix->priorsd,0,priorsd[0]);
   
   for(k=0;k<designmatrix->numparams-1;k++){gsl_vector_set(designmatrix->priormean,k+1,priormean[parentindexes[k]+1]);/** +1 since first entry is constant */
                                            gsl_vector_set(designmatrix->priorsd,k+1,priorsd[parentindexes[k]+1]);/** +1 since first entry is constant */
                                            
   }

/* 
 if(verbose){
  Rprintf("\nCHILD NODE = %s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),nodeid)));
  if(numparents>=1){
        for(k=0;k<numparents;k++){
		Rprintf("PARENT NODE: %s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),parentindexes[k])));}
    Rprintf("\n");		
  } else {Rprintf("No PARENTS\n\n");}
 }
 */
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    dgvalue = gsl_vector_alloc (designmatrix->numparams);/** will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);/** will hold hessian matrix **/
    mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
       
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
 
   /** now we need to solve system defined in laplace_dg()=0 */
  
    iter=0;       
    FDF.f = &laplace_dg;
    FDF.df = &laplace_hessg;
    FDF.fdf = &wrapper_fdf;
    FDF.n = designmatrix->numparams;
    FDF.params = &gparams;
    myBeta = gsl_vector_alloc (designmatrix->numparams);/** this will hold the parameter point estimates */
    /*generate_inits(myBeta,designmatrix);*/     
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    status_inits=generate_inits_n(myBeta,&gparams);
    
    if(status_inits==GSL_SUCCESS){/** model is not singular so find root - optimal **/
  
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
 
   #ifdef PRINTGSL
   print_state (iter, s);
   #endif 
    iter=0; 
       do
         {
           iter++;
       
           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           print_state (iter, s);
           #endif
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
       if( (status != GSL_SUCCESS) && verbose){Rprintf ("Zero finding warning: status = %s at nodeid %d\n", gsl_strerror (status),nodeid);}
       gsl_vector_memcpy(myBeta,s->x);
       
    } /** end of root finding **/     
  /** we now have all the individual parts so put it together to the laplace approx */
  if(status != GSL_SUCCESS){logscore= -DBL_MAX; /** root finding failed so discard model by setting fit to worst possible */
  } else {
  laplace_g(myBeta,&gparams, &gvalue);
  laplace_hessg(myBeta,&gparams, hessgvalue);
   n=obsdata->numDataPts;
   m=designmatrix->numparams;
   perm = gsl_permutation_alloc (m);
   status2=gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
   if(status2 != GSL_SUCCESS){logscore= -DBL_MAX;
   } else {
   logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2)*log((2*M_PI)/n);} /** this is the final value */
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
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp1);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_matrix_free(datamatrix);
   gsl_permutation_free(initsperm);
   if(status == GSL_SUCCESS){gsl_permutation_free(perm);} /** only allocate this is status==GSL_SUCCESS */
   gsl_multiroot_fdfsolver_free (s);
     
   return(logscore); 
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
int laplace_g (const gsl_vector *beta, void *params,double *gvalue)
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
       double term1=0;
       double term2=0;
       double term3=0;
       double storedbl1;
       unsigned int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
     
     /** DO IN THREE PARTS - term1, term2, term3 */
     /** R code "term2<-sum( log(1/(sqrt(2*pi)*sd.loc)) );" **/
     /*term2 = m*(-log(sqrt(2.0*M_PI)*priorsd));*//** assumes priorsd SAME for all parameters => m*() **/
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
int laplace_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues)
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
       double n=Y->size;/** no. observations **/

       unsigned int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
     
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
     
     gsl_vector_memcpy(dgvalues,term1);
    
       
 return GSL_SUCCESS;
     }
       
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues)
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
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** no. params to estimate*/
       double tmp1=0;double tmp2=0;double tmp3=0;

       unsigned int i=0;unsigned int j=0;unsigned int k=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
     
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
                            
                            gsl_blas_ddot (vectmp1long, vectmp2long, gsl_matrix_ptr(hessgvalues,j,k));/** DOT product simply to calcu sum value */
                    } else {*gsl_matrix_ptr(hessgvalues,j,k)=gsl_vector_get(term1,j);}
                     }
                     }
       
 return GSL_SUCCESS;
     }
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_fdf (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_dg(beta, gparams, dgvalues);
       laplace_hessg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/       
void print_state (int iter, gsl_multiroot_fdfsolver * s)
{
 unsigned int i=0;
    Rprintf ("iter = %3u\n",iter);
    
    for(i=0;i< (s->x)->size-1; i++){
          Rprintf ("x=%5.10f ",gsl_vector_get (s->x, i));}
          Rprintf ("x=%5.10f\n",gsl_vector_get (s->x, (s->x)->size-1));
	  
    for(i=0;i< (s->x)->size-1; i++){
          Rprintf ("f(x)=%5.10f ",gsl_vector_get (s->f, i));}
          Rprintf ("f(x)=%5.10f\n",gsl_vector_get (s->f, (s->x)->size-1));   
    
  
}   
   
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
#ifdef OLDinit
int generate_inits(gsl_vector *myBeta,datamatrix *designmatrix){

    /** we want to find the value that maximise the loglikelihood so use
        a very use simple ad-hoc rule: log(p/(1-p)) where p=#success/#trials
        e.g. logit of the proportion of success in response variable **/
    /** get the total number of positive using dot product */

    double cnt=-1.0;
    double p=-1.0;
    double logitP=0.0;
    gsl_vector *Y=designmatrix->Y;
    
    gsl_blas_ddot (Y, Y, &cnt);   /** easy way to do sum(Y) **/
    p=cnt/(Y->size);
    logitP=log(p/(1-p));
    
    gsl_vector_set_all(myBeta,logitP);
   
    #ifdef IND
    Rprintf("manual inits - length myBeta %u\n",myBeta->size);
    gsl_vector_set(myBeta,0,-2.3);
    gsl_vector_set(myBeta,1,0.13);
    gsl_vector_set(myBeta,2,-0.35);
    gsl_vector_set(myBeta,3,0.18);
    gsl_vector_set(myBeta,4,0.14);
    gsl_vector_set(myBeta,5,-0.25);
    gsl_vector_set(myBeta,6,-0.61);
    gsl_vector_set(myBeta,7,-0.17);
    gsl_vector_set(myBeta,8,-0.4);
    gsl_vector_set(myBeta,9,-0.32);
    gsl_vector_set(myBeta,10,-0.22);
    #endif  
    return GSL_SUCCESS;
}
#endif 
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int generate_inits_n(gsl_vector *myBeta,struct fnparams *gparams){

    /** this is the SAME CODE as in the Gaussian case  */
    
    /** beta_hat= (X^T X)^{-1} X^T y **/
    
       const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_permutation *perm = gparams->perm;
     unsigned int i;
     int ss;
     
    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2); */
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */                   
    gsl_linalg_LU_decomp(mattmp3,perm,&ss);
    gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */                   
    
    gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1); /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2); 
    
    for(i=0;i<myBeta->size;i++){gsl_vector_set(myBeta,i,gsl_vector_get(vectmp2,i));} /** set to Least squares estimate */
   
    return GSL_SUCCESS;
}   
  
/** ****************************************************************************************************/   
/** AS ABOVE BUT GAUSSIAN MODEL ************************************************************************/   
/** ****************************************************************************************************
 ***** calc an individual linear Gaussian regression model 
 *******************************************************************************************************/
double calc_Gaussiannode_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix, const double *priormean, const double *priorsd,
                                 const double *priorgamshape, const double *priorgamscale, int gaussiannodeid,const int maxiters, const double epsabs)
{
 int i,ss,status,status2,status_inits;
 unsigned int j,k;
 int iter=0;
 unsigned int numparents=0;
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 
 int *parentindexes=nodescore->parentindexes;/** just memory space **/
 
 gsl_vector *Y,*myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,/**dgvalue,*/*term1,*term2,*term3,*vectmp3long,*vecpriormean,*vecpriorsd,
            *vecpriorgamshape,*vecpriorgamscale,*localbeta;
 gsl_matrix *hessgvalue,/**mattmp1,*/*mattmp2,*mattmp3,*mattmp4,*datamatrix;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 
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
/* 
 if(verbose){
  Rprintf("\nCHILD NODE = %s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),nodeid)));
  if(numparents>=1){
        for(k=0;k<numparents;k++){
		Rprintf("PARENT NODE: %s\n",STRING_PTR(STRING_ELT(VECTOR_ELT(R_labels,0),parentindexes[k])));}
    Rprintf("\n");		
  } else {Rprintf("No PARENTS\n\n");}
 }
 */
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
   /* dgvalue = gsl_vector_alloc (designmatrix->numparams);*//** will hold partial derivates **/
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);/** will hold hessian matrix NOTE: NOT used in solver
                                                                                            as it allocates its own storage based on the problem dimension**/
    /*mattmp1 = gsl_matrix_alloc (obsdata->numparams,obsdata->numparams);*/
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions **/
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
   
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
   /*gparams.mattmp1=mattmp1; */
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.beta=localbeta;
   gparams.perm=initsperm;
 
   /** now we need to solve system defined in laplace_dg()=0 */

    iter=0;       
    FDF.f = &laplace_gaus_dg;
    FDF.df = &laplace_gaus_hessg;
    FDF.fdf = &wrapper_gaus_fdf;
    FDF.n = designmatrix->numparams+1;/** +1 for the precision term */
    FDF.params = &gparams;  
   
    myBeta = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the parameter point estimates INC tau-precision*/
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams+1); /** +1 for the precision term */
    status_inits=generate_gaus_inits(myBeta,&gparams);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    if(status_inits==GSL_SUCCESS){
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
      
    iter=0; 
       do
         {
           iter++;
       
           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           print_state (iter, s);
           #endif
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
     
       if( (status != GSL_SUCCESS) && verbose){Rprintf ("Zero finding warning: status = %s at nodeid %d\n", gsl_strerror (status),nodeid);}
       gsl_vector_memcpy(myBeta,s->x);
            
    } /** end of root finding **/         
  
  /** ************************************/         
  /** we now have all the individual parts so put it together to the laplace approx */
  if(status != GSL_SUCCESS){logscore= -DBL_MAX; /** root finding failed so discard model e.g. set fit to worst possible */
  } else {
  laplace_gaus_g(myBeta,&gparams, &gvalue);
  laplace_gaus_hessg(myBeta,&gparams, hessgvalue);
   n=obsdata->numDataPts;
   m=designmatrix->numparams+1;
   perm = gsl_permutation_alloc (m);
   status2=gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
   if(status2!= GSL_SUCCESS){logscore= -DBL_MAX;
   } else {
   logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2)*log((2*M_PI)/n);} /** this is the final value */
   }
    /*** Last Step before return - free all the gsl vectors, matrices, other etc **/
   gsl_vector_free(Y);
   gsl_vector_free(myBeta);
   gsl_vector_free(localbeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   /*gsl_vector_free(dgvalue);*/
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_vector_free(vecpriormean);
   gsl_vector_free(vecpriorsd);
   gsl_vector_free(vecpriorgamshape);
   gsl_vector_free(vecpriorgamscale);
   gsl_matrix_free(hessgvalue);
   /*gsl_matrix_free(mattmp1);*/
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_matrix_free(datamatrix);
   gsl_permutation_free(initsperm);
   if(status == GSL_SUCCESS){gsl_permutation_free(perm);} /** only allocate this is status==GSL_SUCCESS */
   gsl_multiroot_fdfsolver_free(s);
   
   return(logscore); 
}   
/** GAUSSIAN NODE ************************************************************************************************/  
/** **************************************************************************************************************/
/** g(y) = -(1/n)* log( f(D|betas)f(betas) */ 
/** **************************************************************************************************************/
int laplace_gaus_g (const gsl_vector *betaincTau, void *params,double *gvalue)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
	     /*gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double term2=0;
       double term3=0;
       double term4=0;
       double term5=0;
       double term6=0;
       double storedbl1,storedbl2,storedbl3;
      
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       /*Rprintf("passed tau=%f\n",tau);
       Rprintf("Y=%f\n",gsl_vector_get(Y,0));
       Rprintf("X[0]=%f %f %f\n",gsl_matrix_get(X,0,0),gsl_matrix_get(X,0,1),gsl_matrix_get(X,0,2));
       Rprintf("X[10]=%f %f %f\n",gsl_matrix_get(X,10,0),gsl_matrix_get(X,10,1),gsl_matrix_get(X,10,2));*/
       int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/
       }
     
     /** same as in logistic model */
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
     /*Rprintf("term2 (Rterm3)=%f\n",term2);*/
     
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
     /*Rprintf("term3 (Rterm4)=%f\n",term3); */    
     
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
     /*Rprintf("term6 (Rterm5)=%5.10f\n",term6);*/
     *gvalue=(-1.0/n)*(term2+term3+term4+term5+term6);

     /*Rprintf("loglike=%f\n",term5+term4);*/

       return GSL_SUCCESS;
     }
        
     
/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_dg (const gsl_vector *betaincTau, void *params, gsl_vector *dgvalues)
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
       /*gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1;*/
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       int i=0;
       double storedbl1/*,storedbl2,storedbl3,tmptau*/;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));}
       
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
     for(i=0;i<m;i++){gsl_vector_set(dgvalues,i,gsl_vector_get(term1,i));}
     
     /** dg_dtau **/
     gsl_vector_scale(vectmp1long,-1.0);/** from above vectmp1long is X*beta so this is -X*beta */
     gsl_vector_add(vectmp1long,Y);/** Y-X*beta **/
     gsl_vector_memcpy(vectmp2long,vectmp1long);
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** get sum((Y_i-X_i*beta)^2) using dot product */
     storedbl1= storedbl1*(-0.5);
     gsl_vector_set(dgvalues,m,(-1.0/n)*(
                                          (n/(2.0*tau))+
				          storedbl1 + 
				          (gsl_vector_get(priorgamshape,0)-1.0)/tau
				          - (1.0/gsl_vector_get(priorgamscale,0)))
					  
                                          );      
 return GSL_SUCCESS;
     }


/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_hessg (const gsl_vector *betaincTau, void *params, gsl_matrix *hessgvalues)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
       /*gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*//** numobs long **/
       /*const gsl_vector *priormean = ((struct fnparams *) params)->priormean; */
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       /*const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;*/
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       /*gsl_vector *term2 = ((struct fnparams *) params)->term2;*/
       /*gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
       gsl_matrix *mattmp3 = ((struct fnparams *) params)->mattmp3;
       gsl_matrix *mattmp2 = ((struct fnparams *) params)->mattmp2;
       int n=Y->size;/** no. observations **/
       int m=X->size2;/** number of coefficients excluding tau-precision */
       /*Rprintf("M=%d %f\n",m, gsl_matrix_get(X,2,2));*/
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       int i=0;
       int j,k;
       double tmp1;
      /* double storedbl1,storedbl2,storedbl3,tmptau; */
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));}
 
  /** matrix X%*%t(X) - this is symmetrical so dg_bj_b_k=dg_b_k_b_j **/
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
		      *gsl_matrix_ptr(hessgvalues,j,k)=tmp1;
                    } else {
		      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k)-1.0/(gsl_vector_get(priorsd,j)*gsl_vector_get(priorsd,j)));
		      *gsl_matrix_ptr(hessgvalues,j,k)=tmp1;}
                     }
                     }
   
  
  /** now for dg_dtau second deriv **/
  tmp1=(-1.0/n)*( (-n/(2.0*tau*tau)) - ( (gsl_vector_get(priorgamshape,0)-1)/(tau*tau)) );
  *gsl_matrix_ptr(hessgvalues,m,m)=tmp1;
  /*Rprintf("final cell is %d %d %f\n",m,m,*gsl_matrix_ptr(hessgvalues,m,m));*/
  
  /** now for dg_dtau_dbeta sum(y_i X_ij - X_i beta X_ij) - last ROW of hessian */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     /*Rprintf("beta=%10.15f %10.15f\n",gsl_vector_get(beta,0),gsl_vector_get(beta,1));
     Rprintf("==%f %f\n",gsl_vector_get(vectmp2,0),gsl_vector_get(vectmp2,1));
     Rprintf("==%f %f\n",gsl_vector_get(vectmp1,0),gsl_vector_get(vectmp1,1));*/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_memcpy(term1,vectmp1);/** the remaining term in dg_beta **/
     gsl_vector_scale(term1,(-1.0/n));
     
     /** last row in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvalues,m,j)=gsl_vector_get(term1,j);}
 
 /** dg_dbeta_dtau is same as dg_dtau_dbeta - note symmetry here **/
     
     /** last col in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvalues,j,m)=gsl_vector_get(term1,j);}
 
  return GSL_SUCCESS;
} 
 
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int generate_gaus_inits(gsl_vector *myBeta,struct fnparams *gparams){

    /** just hard coded for moment */
           
    /** use least squares estimates to get started **/
    /** beta_hat= (X^T X)^{-1} X^T y **/
       int status=0;
       const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_vector *vectmp1long = gparams->vectmp1long;/** numobs long **/
       gsl_vector *vectmp2long = gparams->vectmp2long;/** numobs long **/
      /* const gsl_vector *priormean = gparams->priormean;
       const gsl_vector *priorsd   = gparams->priorsd;
       const gsl_vector *priorgamshape = gparams->priorgamshape;
       const gsl_vector *priorgamscale = gparams->priorgamscale; */
       /*gsl_vector *beta = gparams->beta;*/
       /*gsl_vector *term1 = gparams->term1;
       gsl_vector *term2 = gparams->term2;
       gsl_vector *term3 = gparams->term3;*/
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_permutation *perm = gparams->perm;
      double n=Y->size;/** no. observations **/
     double m=X->size2;/** number of coefficients excluding tau-precision */
     unsigned int i;
     int ss;
     double variance=0.0;
     
    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2); */
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */                   
    status=gsl_linalg_LU_decomp(mattmp3,perm,&ss);
    if(status!=GSL_SUCCESS){return(GSL_FAILURE);}
    status=gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */                   
    if(status!=GSL_SUCCESS){return(GSL_FAILURE);}
    gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1); /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2); 
    
    for(i=0;i<myBeta->size-1;i++){gsl_vector_set(myBeta,i,gsl_vector_get(vectmp2,i));} /** set to Least squares estimate */
   
    /** now for variance estimate */
    /** first get y_hat estimate */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp2, 0.0, vectmp1long); /** vectmp1 is y_hat */ 
    gsl_vector_scale(vectmp1long,-1.0);/** - y_hat */
    gsl_vector_add(vectmp1long,Y);/** now have Y-y_hat (or -y_hat + Y) */
    gsl_vector_memcpy(vectmp2long,vectmp1long);
    gsl_blas_ddot (vectmp1long, vectmp2long, &variance);/** got sum((Y-Y_hat)^2) */
    variance=variance/(n-m);/** unbiased estimator using denominator n-#term in regression equation **/
    
    gsl_vector_set(myBeta,myBeta->size-1,1.0/variance);/** as are using precision parameterization not variance **/ 
    
   /* Rprintf("b0=%f b1=%f b2=%f\n",gsl_vector_get(myBeta,0),gsl_vector_get(myBeta,1),gsl_vector_get(myBeta,2));
    Rprintf("precision=%f\n",gsl_vector_get(myBeta,3));*/
   /* double mean=gsl_stats_mean(Y,1,Y->size);
    double var=gsl_stats_variance(Y,1,Y->size);
    int i=0;
    gsl_vector_set(myBeta,0,mean);
    for(i=1;i<myBeta->size-1;i++){gsl_vector_set(myBeta,i,mean);}
    gsl_vector_set(myBeta,myBeta->size-1,0.5);
     */
    
    return GSL_SUCCESS;
}   
 

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_gaus_fdf (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_gaus_dg(beta, gparams, dgvalues);
       laplace_gaus_hessg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }   
     

 
   
   
   
   
   
   
   
   
   
   
   
   
   
