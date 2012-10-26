/** **********************************************************************************************************************/
/** **********************************************************************************************************************/
/** README. */
/** The top level function in this file is calc_node_Score_binary_rv(), which computes the log marginal likelihood for a */
/** binary logit model with one level of random effects to adjust for within group correlation. Due to the numerics are  */
/** far more involved than for glm's - in particular the Laplace approx to the mlik involves summing over Laplace approxs*/
/** - one for each group of observations. Hence there is a single outer Laplace calc. and then many "inner" Laplace calcs*/
/** key differences with non-mixed models - the design matrix contains the local rv term (epsilon, say). There are a     */ 
/** number of different accuracy parameters to deal with e.g. step sizes for numerical derivs etc.                       */
/** init values - as for binary use glm least square ests but do NOT use Gaussian variance value but rather the use      */
/** p(1-p) where p is the intercept term in the mean e.g. p=)(exp(b0)/(1+exp(b0))                                        */

/** top level - calc_node_Score_binary_rv(): uses a minimiser (bfgs2) rather than root finder and 
                                  g_outer(): the objective function -1/n*log(f(D|theta)f(thera))
                              rv_dg_outer(): the derivative of the objective function - note this is a numeric derivative
                                             using central finite differences
                           rv_hessg_outer(): numerical derivative using central finite differences
                           
    second level -                g_outer(): this is basically just a wrapper which computes the prior for the means and
                                             calls a function which computes a Laplace approx for each separate data group 
                                  g_inner(): this does the actual Laplace approx calc for each group and uses root 
                                             finding with analytical hessian
    
    third level                   g_inner(): does root finding using next two functions and then computes the Laplace value 
                              rv_dg_inner(): first derviative 
                           rv_hessg_inner(): second derivatives
                               rv_g_inner(): value of the -1/n*log(f(D|theta)f(theta) for each group of data D            
    
    others
                           g_outer_single(): as g_outer() but where one variable is allowed to vary - used in num derivative
                           
                                                                                                                         */
/** **********************************************************************************************************************/
/** **********************************************************************************************************************/
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h" 
#include "node_binomial_rv_inner.h"
#include "node_binomial_rv_Rsupp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#define NEWHESS
#define NOPRIOR1 /**turn on for removing effect of priors - this is purely to enable a check of loglike values with lme4 */
/** ****************************************************************************************************
 ***** calc an individual logistic regression model 
 *******************************************************************************************************/
void calc_node_Score_binary_rv_R(network *dag, datamatrix *obsdata, int nodeid,  int errverbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,int storeModes, double epsabs_inner, int maxiters_inner, double finitestepsize, int verbose,
				 double h_guess, double h_epsabs, int maxiters_hessian)
{
#ifdef NOPRIOR
Rprintf("############ Warning - Priors turned off - use only for checking mlik value! ################\n");
#endif
  
  int i,status,actual_status,sss,index=0,iter;
  /*int j;*/
  gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*localbeta,*localbeta2,*dgvalues,*finitefactors,*factorindexes, *finitestepsize_vec,*nmstepsize;/** this will hold the parameter point estimates + 1 is for precision of rv's */
  struct fnparams gparams;/** for passing to the gsl zero finding functions */
  double gvalue;double nm_size=0.0;
  gsl_matrix *mattmp2,*mattmp3,*mattmp4,*hessgvalues,*hessgvalues3pt;
  double mydet=0.0,logscore=0.0;/*,logscore3pt=0.0;*/
  gsl_permutation *initsperm;
  gsl_permutation *perm=0; 
  const gsl_multimin_fminimizer_type *T;
       gsl_multimin_fminimizer *s;
     gsl_multimin_function F;  
  int n,m;
 /* double min_error,cur_error,accurate_logscore=0,accurate_logscore3pt=0,bestsize=0,lowerend,upperend,h_guess,h_epsabs;*/
  /*const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;*/ 
  /*double h_lowerbound[1],h_upperbound[1],h_guess_array[1];
  int h_nbd[1];*/
  int nDim;/** dimension of optim problem */
  int *nbd;/** nbd is an integer array of dimension nDim.
	                                      On entry nbd represents the type of bounds imposed on the variables, and must be specified as follows:
	                                      nbd(i)=0 if x(i) is unbounded,
		                              1 if x(i) has only a lower bound,
		                              2 if x(i) has both lower and upper bounds, and
		                              3 if x(i) has only an upper bound.
	                                      On exit nbd is unchanged.*/
  
  double *lowerbounds,*upperbounds; /* h_gvalue;*//*,lowestHesserror,beststepsize;*/
  int failcode;/** check code see R ?optim - if non-zero the a problem **/
  double factr=1e-07;/** error size scaler - this is the default value*/
  double pgtol=1e-07;/** default value is zero - this is the gradient tolerance - mmm what does that actually mean? */
  int fncount,grcount;/** hold number of evaluations */
  char msg[60];/** error message */
  int trace=errverbose;/** like verbose */
  int nREPORT=1000;/** report freq*/
  int lmm=5;/** see R ?optim - number of function evals to store - default is 5 */
  /** want to find the modes of the function g(betas) where betas=b_0,b_1,,...,tau, the latter being precision */
  /** g(betas) is not differentiable as it contains integrals (which themselves need Laplace estimates **/
  
    
  /** SETUP things which are the same across all data chunks - groups  */
  /** build design matrix which is designmatrix->datamatrix=X, designmatrix->Y=Y plus priors designmatrix->priorsd, designmatrix->priormean **/
  /** NOTE: design matrix here does include the random effect term **/
  /** note - numparams does NOT include precision term - numpars +1 */
  build_designmatrix_rv(dag,obsdata,priormean, priorsd,priorgamshape,priorgamscale,designmatrix,nodeid,storeModes);
  
  nDim=designmatrix->numparams+1; 
  lowerbounds=(double *)R_alloc(nDim,sizeof(double*));
  upperbounds=(double *)R_alloc(nDim,sizeof(double*));
  nbd=(int *)R_alloc(nDim,sizeof(int*));
  for(i=0;i<nDim-1;i++){lowerbounds[i]=-DBL_MAX;
                        upperbounds[i]=DBL_MAX;
			nbd[i]=0;}
			nbd[nDim-1]=1;lowerbounds[nDim-1]=0.001;/** lower bound for precision */
  finitefactors = gsl_vector_alloc(7);/** used to change stepsize in hessian estimate **/			
  gsl_vector_set(finitefactors,0,1.0E-03);gsl_vector_set(finitefactors,1,1.0E-02);gsl_vector_set(finitefactors,2,1.0E-01);
  gsl_vector_set(finitefactors,3,1.0);gsl_vector_set(finitefactors,4,1.0E+01);gsl_vector_set(finitefactors,5,1.0E+02);
  gsl_vector_set(finitefactors,6,1.0E+03);
  
  factorindexes = gsl_vector_alloc(7);/** used to change stepsize in hessian estimate **/			
  for(i=0;i<7;i++){gsl_vector_set(factorindexes,i,i);}
  
  /** change finite.step.size by 0.1,1, and 10 factors respectively **/
  
  dgvalues = gsl_vector_alloc (designmatrix->numparams+1);/** inc rv precision */
  
  myBeta = gsl_vector_alloc (designmatrix->numparams+1);/** inc rv precision */
  vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
  vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
  mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
  mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
  mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
  initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
  vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
  vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
  localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions - excl. precision **/
  localbeta2 = gsl_vector_alloc (designmatrix->numparams+1);/** scratch space in later functions - inc. precision **/
  
  hessgvalues = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);
  hessgvalues3pt = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);
  
  gparams.designdata=designmatrix;
  
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.perm=initsperm;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.beta=localbeta;
   gparams.betaincTau=localbeta2;
   gparams.epsabs_inner=epsabs_inner;
   gparams.maxiters_inner=maxiters_inner;
   gparams.verbose=verbose;
   gparams.finitestepsize=finitestepsize;
   
   dag->nodeScoresErrCode[nodeid]=0;/** reset error code to no error **/
   
   /*status=GSL_SUCCESS;*/
   generate_rv_inits(myBeta,&gparams);
   /*Rprintf("starting optimisation\n");*/
   /** run a loop over different stepsize - starting with the smallest first as this is more likely successful **/
   for(i=0;i<finitefactors->size;i++){/** go forwards through the factors so start with SMALLEST STEPSIZE */
   /*Rprintf("step size iteration %d\n",i);*/
     failcode=0;/** reset*/
    gparams.finitestepsize=gsl_vector_get(finitefactors,i)*finitestepsize;
   
     lbfgsb(nDim, lmm, myBeta->data, lowerbounds, upperbounds, nbd, &gvalue, &g_outer_R,
                      &rv_dg_outer_R, &failcode, 
	              &gparams,
	              factr,
                      pgtol, &fncount, &grcount,
                      maxiters, msg, trace, nREPORT);
		      
    if(!failcode){dag->nodeScoresErrCode[nodeid]=0;/*bestsize=gparams.finitestepsize;*/break;}/** break out of for loop if no error as we are done **/	     
   
   } /** end of for loop so now have mode estimates */
     
   if(failcode){Rprintf("%s at node %d\n",msg,nodeid+1);/** notify if there is an error and set final error code **/
		     dag->nodeScoresErrCode[nodeid]=1;
   } 
     
    gparams.finitestepsize=finitestepsize;/** reset */
    if(storeModes){/** keep a copy of the parameter modes found for use later in other function calls etc**/
	 index=0;    /*Rprintf("size of beta=%d %f %f\n",myBeta->size, gsl_vector_get(myBeta,0),gsl_vector_get(myBeta,1));*/
		     for(i=0;i<dag->numNodes+3;i++){/** roll myBeta into dag->modes into the appropriate columns**/
		       if(gsl_matrix_get(dag->modes,nodeid,i)!=DBL_MAX){
			 gsl_matrix_set(dag->modes,nodeid,i,gsl_vector_get(myBeta,index++));}} 
                   /*for(i=0;i<dag->numNodes+3;i++){Rprintf("%e ",gsl_matrix_get(dag->modes,nodeid,i));}Rprintf("\n");*/
		   
		   }     

   /** now compute the hessian at the step size with lowest error **/
   /*Rprintf("starting hessian estimation\n");*/
   n=obsdata->numDataPts;
   m=designmatrix->numparams+1;/** inc precision */
   perm = gsl_permutation_alloc (m);
 
   /** just re-use as much of existing gparams as possible - so names are meaningless e.g. betafixed is actually gvalue */
   gparams.betaincTau=myBeta;
   gparams.nDim=n;
   gparams.mDim=m;
   gparams.perm=perm;
   gparams.mattmp2=hessgvalues;
   gparams.mattmp3=hessgvalues3pt;
   gparams.betafixed=gvalue;
   
   
    F.f = &compute_mlik_nm;
    F.params = &gparams;
    F.n = 1;
   
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);
   
    finitestepsize_vec = gsl_vector_alloc (1);
    gsl_vector_set (finitestepsize_vec, 0, h_guess);
    nmstepsize = gsl_vector_alloc (1);
    gsl_vector_set_all (nmstepsize, h_guess); 
    gsl_multimin_fminimizer_set (s, &F, finitestepsize_vec, nmstepsize);
    status = GSL_SUCCESS;
    
     iter=0;
   
    do
         {
           iter++;
           status = gsl_multimin_fminimizer_iterate (s);
     
           if (status) 
             break;
	   
	   nm_size = gsl_multimin_fminimizer_size (s);
           status = gsl_multimin_test_size (nm_size, h_epsabs);
     /*
           if (status == GSL_SUCCESS)
             {
               Rprintf ("converged to minimum at\n");
             }
     */
           Rprintf ("iter=%5d error in mlik=%3.5e using fin.diff step= %3.2e nmsize=%3.2e\n", iter,s->fval,gsl_vector_get (s->x, 0),nm_size);
    
         }
       while (status == GSL_CONTINUE && iter < maxiters_hessian);
       if( (status != GSL_SUCCESS)){actual_status=status;/** copy for use later **/
                                    status=GSL_FAILURE;} /** solution failed to achieve a value below h_epsabs **/                                                               
	 
    finitestepsize=gsl_vector_get(s->x,0);/** get best fin.diff stepsize **/
    
    dag->hessianError[nodeid]= s->fval;/** get fin.diff error **/
    
    gsl_multimin_fminimizer_free (s);
     
   
      
      /*Rprintf("GVALUE=%e nodeid=%d\n",gvalue,nodeid+1);*/
       /*finitestepsize=1E-04;
       h_guess=finitestepsize;
*/
       switch(status){  /** choose which type of node we have */
                     case GSL_SUCCESS:{    
		                     /** successful finite diff so now do final computation with the optimal step size **/
                                     /*Rprintf("search for optimal step size : status = %s at nodeid %d\n", gsl_strerror (status),nodeid+1);*/
                                     rv_hessg_outer(myBeta,&gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** EDIT BACK to "finitestepsize" start with LARGEST STEPSIZE **/
				    /* Rprintf("HESSIAN using stepsize=%e\n",finitestepsize);
				     for(i1=0;i1<hessgvalues->size1;i1++){
				        for(i2=0;i2<hessgvalues->size2;i2++){Rprintf("%e ",gsl_matrix_get(hessgvalues,i1,i2));}Rprintf("\n");}   */
                                     status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
                                     mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
                                     logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
				     if(gsl_isnan(logscore)){logscore= R_NaN;dag->nodeScoresErrCode[nodeid]=2;}
				     dag->nodeScores[nodeid]=logscore;
				       
		                      break;  
		     }
       
		     case GSL_FAILURE: {/** the minimiser did not find a minimum meeting the accuracy requirements and so may be unreliable **/
		                       Rprintf ("-- ERROR! -- search for optimal step size error: status = %s at nodeid %d\n", gsl_strerror (actual_status),nodeid+1);
                                       rv_hessg_outer(myBeta,&gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
                                       /*Rprintf("HESSIAN using stepsize=%e\n",finitestepsize);
				       for(i1=0;i1<hessgvalues->size1;i1++){
				        for(i2=0;i2<hessgvalues->size2;i2++){Rprintf("%e ",gsl_matrix_get(hessgvalues,i1,i2));}Rprintf("\n");} */
				        
				        status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
                                       mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
                                       logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
				       dag->nodeScoresErrCode[nodeid]=4;
				       if(gsl_isnan(logscore)){logscore= R_NaN;dag->nodeScoresErrCode[nodeid]=2;}
				       dag->nodeScores[nodeid]=logscore;
				       
		                       
				       break; 
		     }
		     
		     default:{Rprintf("got case %s\n",gsl_strerror (status)); error("in default switch in calc_node_Score_binary_rv_R() - should never get here!");}  
		     
          }
          
        
   /** try the bounded search for h stepsize rather than one-dim min which needs bound specified **/     
        
  
   /** now free up allocated memory **/
   for(i=0;i<designmatrix->numUnqGrps;i++){gsl_matrix_free(designmatrix->array_of_designs[i]);
                                           gsl_vector_free(designmatrix->array_of_Y[i]);}
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);
   gsl_vector_free(designmatrix->priorgamshape);
   gsl_vector_free(designmatrix->priorgamscale);
   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix_noRV);
   gsl_vector_free(dgvalues);
   gsl_vector_free(myBeta); 
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_permutation_free(initsperm);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(localbeta);
   gsl_vector_free(localbeta2);
   gsl_matrix_free(hessgvalues);
   gsl_matrix_free(hessgvalues3pt);
   gsl_vector_free(finitefactors);
   gsl_vector_free(factorindexes);
   
   /*if(!failcode){*/gsl_permutation_free(perm);/*}*/
   
   /*dag->nodeScores[nodeid]=logscore;*/

}
/** **************************************************************************************************************/
/** **************************************************************************************************************/
 
/** *************************************************************************************************************/
/** **************************************************************************************************************/
/** g_outer(y) = -(1/n)* log( prod[g_inner()].prior )                                                            */ 
/** this function is the product (or sum of logs) of the likelihoods for each individual groups of data          */
/** where these individual likelihood have integrals in them - to integrat out the latent rv effect. g_outer     */
/** then calls g_inner() which then deals with each of the individual likelihoods in turn                        */
/** **************************************************************************************************************/
/** **************************************************************************************************************/
double g_outer_R (int Rn, double *betaincTauDBL, void *params) /*typedef double optimfn(int n, double *par, void *ex);*/
{
  int i,j;
  double term1=0.0,singlegrp=0.0;
  const datamatrix *designdata = ((struct fnparams *) params)->designdata;/** all design data inc Y and priors **/

const gsl_vector *priormean = designdata->priormean;
  const gsl_vector *priorsd   = designdata->priorsd;
  const gsl_vector *priorgamshape   = designdata->priorgamshape;
  const gsl_vector *priorgamscale   = designdata->priorgamscale;
   gsl_vector *beta   = ((struct fnparams *) params)->beta;/** does not include precision */
  gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
  gsl_vector *vectmp2 =((struct fnparams *) params)->vectmp2;/** numparams long*/
  gsl_vector *betaincTau=((struct fnparams *) params)->betaincTau;/** to copy betaincTauDBL into **/
  double epsabs_inner=((struct fnparams *) params)->epsabs_inner;/** absolute error in internal laplace est */
  int maxiters_inner=((struct fnparams *) params)->maxiters_inner;/** number of steps for inner root finder */
  int verbose=((struct fnparams *) params)->verbose;/**  */
  
  int n_betas= (designdata->datamatrix_noRV)->size2;/** number of mean terms excl rv and precision **/
  int n=(designdata->datamatrix_noRV)->size1;/** total number of obs **/
  
  double term2=0.0,term3=0.0,term4=0.0,gval=0.0;
  /*Rprintf("%d %d\n",n_betas,betaincTau->size);*/
  double tau;
  
  for(i=0;i<betaincTau->size;i++){gsl_vector_set(betaincTau,i,betaincTauDBL[i]);} /** copy R double array into gsl vect **/
  
  /*Rprintf("got = %f %f %f\n",gsl_vector_get(betaincTau,0),gsl_vector_get(betaincTau,1),gsl_vector_get(betaincTau,2));*/
    
  tau=gsl_vector_get(betaincTau,n_betas);/** extract the tau-precision from *beta - last entry */
  /*Rprintf("g_outer_ tau=%f\n",tau);*/
 
  if(tau<0.0){Rprintf("tau negative in g_outer!\n");error("");}
  
  /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<n_betas;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/
       }
     
  /** part 1 - the integrals over each group of observations - use laplace for this and that is taken care of in g_inner */ 
  /** first we want to evaluate each of the integrals for each data group **/ 
       for(j=0;j<designdata->numUnqGrps;j++){/** for each data group **/
	/*j=0;*/
	/* Rprintf("processing group %d\n",j+1);
	 Rprintf("tau in loop=%f\n",gsl_vector_get(betaincTau,n_betas));*/
	  singlegrp=g_inner(betaincTau,designdata,j, epsabs_inner,maxiters_inner,verbose);
	 
	 if(gsl_isnan(singlegrp)){error("nan in g_inner\n");}
	  term1+= singlegrp;
      }
      
/*Rprintf("term1 in g_outer=%f\n",term1);*/	
  /** part 2 the priors for the means **/
  term2=0; for(i=0;i<n_betas;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
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
  
  
  /** part 3 the prior for the precision tau **/
  term4=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau)
	     -(tau/gsl_vector_get(priorgamscale,0));
   
	     
   gval=(-1.0/n)*(term1+term2+term3+term4);
   /** NO PRIOR */
  /* Rprintf("WARNING - NO PRIOR\n");*/
  #ifdef NOPRIOR
  gval=(-1.0/n)*(term1);
  #endif
   if(gsl_isnan(gval)){error("g_outer_R\n");}
/*Rprintf("g_outer_final=%f term1=%f term2=%f term3=%f term4=%f total=%f %d\n",gval,term1,term2,term3,term4,term1+term2+term3+term4,n);	*/
	return(gval);/** negative since its a minimiser */
}	
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/ 
void rv_dg_outer_R (int n, double *betaDBL, double *dgvaluesDBL,void *params)
{
  struct fnparams *gparams = ((struct fnparams *) params);
  
  /** fixed covariate and precision terms **/
  
  gsl_function F;int i;int haveTau=0;
  double result, abserr;
  double h=((struct fnparams *) gparams)->finitestepsize;
  /*double h_adj=0.0;*//** a new h is existing h is too small **/
  gsl_vector *betaincTau=((struct fnparams *) gparams)->betaincTau;/** scratch space to copy betaincTauDBL into **/
  
  for(i=0;i<betaincTau->size;i++){gsl_vector_set(betaincTau,i,betaDBL[i]);} /** copy into gsl_vector **/
   
  if(betaDBL[n-1]<0.0){error("negative tau in rv_dg_outer_R\n");}
  /** if get tau which is negative */
  
       F.function = &g_outer_single;
       F.params = gparams;
   
  for(i=0;i<n;i++){ /** each of the partial derivatives */    
  gparams->fixed_index=i;
  if(i== n-1){haveTau=1;} else {haveTau=0;}
    
  /** readme - evaluating f() at a negative value e.g. tau-h **/
  if(!haveTau){gsl_deriv_central (&F, betaDBL[i], h, &result, &abserr);/*Rprintf("fixed=%d val=%f\n",i,result);*/
  } else { /** first try central and if this goes into negative tau area then revert to forwards **/
           gsl_deriv_central(&F, betaDBL[i], h, &result, &abserr); 
	   if(gsl_isnan(abserr)){gsl_deriv_forward(&F, betaDBL[i], h, &result, &abserr);}
  }
  
  dgvaluesDBL[i]=result;
  }
  
 
  
}
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
int rv_hessg_outer( gsl_vector* beta, void* params, gsl_matrix* hessgvalues,double h, gsl_matrix* hessgvalues3pt)
{
  
  struct fnparams *gparams = ((struct fnparams *) params);
  int i,j;
  int haveTau=0;
  gsl_function F;
  gparams->betaincTau=beta;  
  
   
  /** if get tau which is negative */
  F.function = &g_outer_single;
  F.params = gparams;
  
 if(gsl_vector_get(beta,beta->size-1)<0.0){Rprintf("negative tau in hess %e\n",gsl_vector_get(beta,beta->size-1));error("");}
 
  /** diagnonal terms d^2f/dx^2 - finite differences - difference between two */
  for(i=0;i<hessgvalues->size1;i++){
   for(j=0;j<hessgvalues->size2;j++){
        if(i>=j){/** hessian is symmetric **/
	        /** d^2f/ dbeta[i] dbeta[j] */
	 gparams->fixed_index=i;/** d^2f/dbeta1 dbeta0 = central diff of beta1 when beta0=beta0+h/2 and similar when beta0=beta0-h/2  -rest are fixed */
         if(i== hessgvalues->size1-1){haveTau=1;} else {haveTau=0;}
   
         
	 gsl_matrix_set(hessgvalues,i,j,get_second_deriv_5pt(gparams,i,j,h,haveTau,&F));/** last args: havetau=0 **/
	 gsl_matrix_set(hessgvalues3pt,i,j,get_second_deriv_3pt(gparams,i,j,h,haveTau,&F));/** last args: havetau=0 **/
	} /** end i>j **/
	
   }
  }
  
  /** finally fill out matrix on upper diagnonal */
  for(i=0;i<hessgvalues->size1;i++){
   for(j=0;j<hessgvalues->size2;j++){
     if(i>=j){gsl_matrix_set(hessgvalues,j,i,gsl_matrix_get(hessgvalues,i,j));}}}
 
 /** finally fill out matrix on upper diagnonal */
  for(i=0;i<hessgvalues3pt->size1;i++){
   for(j=0;j<hessgvalues3pt->size2;j++){
     if(i>=j){gsl_matrix_set(hessgvalues3pt,j,i,gsl_matrix_get(hessgvalues3pt,i,j));}}}
 
  return GSL_SUCCESS;
  
}

/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double get_second_deriv_5pt(struct fnparams *gparams, int i, int j, double h, int haveTau, gsl_function *F)
{
  double result1,result2,result3, result4, result5, abserr1,abserr2,abserr3,abserr4, abserr5;
  
  gsl_vector *beta=gparams->betaincTau;
  double *beta_j=&(beta->data[j]);/** pointers to the relevant enties in beta vector **/
  double *beta_i=&(beta->data[i]);
  const double masterbetaj=gsl_vector_get(beta,j);
  /** want to call g_outer_single with varible j shifted */
  
  if(!haveTau){/** if not tau  use central differences **/
    
    /** 4 terms for df_xj each of which need expanded. We FIX x_j at different value and then expand five point formula on xi **/
    
    /** f(x_j-2h,x_i...) etc */
    *beta_j=*beta_j-2.0*h;
    gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);/** this is value of first derivative at b_j=b_j-2h */
    *beta_j=masterbetaj;/** reset **/
    
     /** f(x_j-h,x_i...) etc */
    *beta_j=*beta_j-1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
    *beta_j=masterbetaj;/** reset **/
    
    /** f(x_j+h,x_i...) etc */
    *beta_j=*beta_j+1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result3, &abserr3);
    *beta_j=masterbetaj;/** reset **/	
    
    /** f(x_j+2h,x_i...) etc */
    *beta_j=*beta_j+2.0*h;
    gsl_deriv_central(F, *beta_i, h, &result4, &abserr4);	   
    *beta_j=masterbetaj;/** reset **/
  
  return((1.0/(12.0*h))*(result1-8.0*result2+8.0*result3-result4));
  }
  
  if(haveTau && i==j && *beta_i-2.0*h<0.0){/** want d^2f/dtau dtau and tau would be negative given using a central diff so use left end version */
    
     /** f(x_j,x_i...) etc */
     *beta_j=*beta_j;/** no change */
     gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);
     if(gsl_isnan(abserr1)){gsl_deriv_forward (F, *beta_i, h, &result1, &abserr1);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/  
     
      /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
     if(gsl_isnan(abserr2)){gsl_deriv_forward (F, *beta_i, h, &result2, &abserr2);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/
     
     /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+2.0*h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result3, &abserr3);
     if(gsl_isnan(abserr3)){gsl_deriv_forward (F, *beta_i, h, &result3, &abserr3);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/
     
     /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+3.0*h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result4, &abserr4);
     if(gsl_isnan(abserr4)){gsl_deriv_forward (F, *beta_i, h, &result4, &abserr4);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/ 
     
     /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+4.0*h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result5, &abserr5);
     if(gsl_isnan(abserr5)){gsl_deriv_forward (F, *beta_i, h, &result5, &abserr5);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/ 
     
     return((1.0/(12.0*h))*(-25.0*result1+48.0*result2-36.0*result3+16.0*result4-3.0*result5));
  }
  
  if(haveTau){/** want d^2f/dtau dx or  d^2f/dtau dtau and in the latter we can evalute use a central difference for the first derivative */
  
    /** f(x_j-2h,x_i...) etc */
    *beta_j=*beta_j-2.0*h;
    gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);/** this is value of first derivative at b_j=b_j-2h */
    if(gsl_isnan(abserr1)){gsl_deriv_forward (F, *beta_i, h, &result1, &abserr1);}
    *beta_j=masterbetaj;/** reset **/
    
     /** f(x_j-h,x_i...) etc */
    *beta_j=*beta_j-1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
    if(gsl_isnan(abserr2)){gsl_deriv_forward (F, *beta_i, h, &result2, &abserr2);}
    *beta_j=masterbetaj;/** reset **/
    
    /** f(x_j+h,x_i...) etc */
    *beta_j=*beta_j+1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result3, &abserr3);
    if(gsl_isnan(abserr3)){gsl_deriv_forward (F, *beta_i, h, &result3, &abserr3);}
    *beta_j=masterbetaj;/** reset **/	
    
    /** f(x_j+2h,x_i...) etc */
    *beta_j=*beta_j+2.0*h;
    gsl_deriv_central(F, *beta_i, h, &result4, &abserr4);
    if(gsl_isnan(abserr4)){gsl_deriv_forward (F, *beta_i, h, &result4, &abserr4);}	   
    *beta_j=masterbetaj;/** reset **/

  return((1.0/(12.0*h))*(result1-8.0*result2+8.0*result3-result4));
  }
  
  error("should never get here - hessian\n");
  return(1.0);
}
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double get_second_deriv_3pt(struct fnparams *gparams, int i, int j, double h, int haveTau, gsl_function *F)
{
  double result1,result2,result3, abserr1,abserr2,abserr3;
  
  gsl_vector *beta=gparams->betaincTau;
  double *beta_j=&(beta->data[j]);/** pointers to the relevant enties in beta vector **/
  double *beta_i=&(beta->data[i]);
  const double masterbetaj=gsl_vector_get(beta,j);
  /** want to call g_outer_single with varible j shifted */
  
  if(!haveTau){/** if not tau  use central differences **/
    
    /** 2 terms for df_xj each of which need expanded. We FIX x_j at different value and then expand five point formula on xi **/
    
     /** f(x_j-h,x_i...) etc */
    *beta_j=*beta_j+1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);
    *beta_j=masterbetaj;/** reset **/
    
    /** f(x_j+h,x_i...) etc */
    *beta_j=*beta_j-1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
    *beta_j=masterbetaj;/** reset **/	
  
  return((1.0/(2.0*h))*(result1-result2));
  }
  
  if(haveTau && i==j && *beta_i-1.0*h<0.0){/** want d^2f/dtau dtau and tau would be negative given using a central diff so use left end version */
    
     /** f(x_j,x_i...) etc */
     *beta_j=*beta_j;/** no change */
     gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);
     if(gsl_isnan(abserr1)){gsl_deriv_forward (F, *beta_i, h, &result1, &abserr1);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/  
     
      /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
     if(gsl_isnan(abserr2)){gsl_deriv_forward (F, *beta_i, h, &result2, &abserr2);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/
     
     /** f(x_j+h,x_i...) etc */
     *beta_j=*beta_j+2.0*h;/** no change */  
     gsl_deriv_central(F, *beta_i, h, &result3, &abserr3);
     if(gsl_isnan(abserr3)){gsl_deriv_forward (F, *beta_i, h, &result3, &abserr3);} /** in case h used trips into negative tau */
     *beta_j=masterbetaj;/** reset **/
     
     return((1.0/(2.0*h))*(-3.0*result1+4.0*result2-result3));
  }
  
  if(haveTau){/** want d^2f/dtau dx or  d^2f/dtau dtau and in the latter we can evalute use a central difference for the first derivative */
  
    /** f(x_j-2h,x_i...) etc */
    *beta_j=*beta_j+1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result1, &abserr1);/** this is value of first derivative at b_j=b_j-2h */
    if(gsl_isnan(abserr1)){gsl_deriv_forward (F, *beta_i, h, &result1, &abserr1);}
    *beta_j=masterbetaj;/** reset **/
    
     /** f(x_j-h,x_i...) etc */
    *beta_j=*beta_j-1.0*h;
    gsl_deriv_central(F, *beta_i, h, &result2, &abserr2);
    if(gsl_isnan(abserr2)){gsl_deriv_forward (F, *beta_i, h, &result2, &abserr2);}
    *beta_j=masterbetaj;/** reset **/
  

  return((1.0/(2.0*h))*(result1-result2));
  }
  
  error("should never get here - hessian\n");
  return(1.0);
}
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double compute_mlik(double finitestepsize, void *params)
{
   struct fnparams *gparams = ((struct fnparams *) params);
   gsl_vector *myBeta=gparams->betaincTau;
   int n=gparams->nDim;
   int m=gparams->mDim;
   gsl_permutation *perm=gparams->perm;
   gsl_matrix *hessgvalues=gparams->mattmp2;
   gsl_matrix *hessgvalues3pt=gparams->mattmp3;
   double gvalue=gparams->betafixed;
   int status,sss;
   double mydet;
   double logscore,logscore3pt;
   /*int i,j;*/
   
   /*Rprintf("got h=%e n=%d m=%d gvalue=%e\n",finitestepsize,n,m,gvalue);
  for(i=0;i<myBeta->size;i++){Rprintf("beta= %f ",gsl_vector_get(myBeta,i));}Rprintf("\n");
  */ 
   
   rv_hessg_outer(myBeta,gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
   
   /*for(i=0;i<hessgvalues3pt->size1;i++){for(j=0;j<hessgvalues3pt->size2;j++){Rprintf("%e ",gsl_matrix_get(hessgvalues3pt,i,j));}Rprintf("\n");}*/
   
   status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
  
   status=gsl_linalg_LU_decomp(hessgvalues3pt,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues3pt);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore3pt= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
   
   /*gparams->logscore=logscore;
   gparams->logscore3pt=logscore3pt;*/
   /*Rprintf("logscore=%e logscore3pt=%e\n",logscore,logscore3pt);*/
   return(fabs(logscore-logscore3pt)); 
 
}

/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double compute_mlik_nm(const gsl_vector *finitestepsize_vec, void *params)
{
   struct fnparams *gparams = ((struct fnparams *) params);
   gsl_vector *myBeta=gparams->betaincTau;
   int n=gparams->nDim;
   int m=gparams->mDim;
   gsl_permutation *perm=gparams->perm;
   gsl_matrix *hessgvalues=gparams->mattmp2;
   gsl_matrix *hessgvalues3pt=gparams->mattmp3;
   double gvalue=gparams->betafixed;
   int status,sss;
   double mydet;
   double logscore,logscore3pt;
   double finitestepsize=gsl_vector_get(finitestepsize_vec,0);
   
   /*int i,j;
   
   Rprintf("got h=%e n=%d m=%d gvalue=%e\n",finitestepsize,n,m,gvalue);
  for(i=0;i<myBeta->size;i++){Rprintf("beta= %f ",gsl_vector_get(myBeta,i));}Rprintf("\n");
   */
   
   rv_hessg_outer(myBeta,gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
   
   /*for(i=0;i<hessgvalues3pt->size1;i++){for(j=0;j<hessgvalues3pt->size2;j++){Rprintf("%e ",gsl_matrix_get(hessgvalues3pt,i,j));}Rprintf("\n");}*/
   
   status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
  
   status=gsl_linalg_LU_decomp(hessgvalues3pt,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues3pt);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore3pt= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
   
   /*gparams->logscore=logscore;
   gparams->logscore3pt=logscore3pt;*/
   /*Rprintf("logscore=%e logscore3pt=%e\n",logscore,logscore3pt);*/
   return(fabs(logscore-logscore3pt)); 
 
}
