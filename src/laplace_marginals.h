#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
     
void calc_network_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix,  const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale,
				 int nodenum, int varnum, int whichgaus,gsl_matrix *posterior,const int maxiters, const double epsabs);

void calc_node_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
             datamatrix *designmatrix,  const double *priormean, const double *priorsd, const double *priorgamshape,
             const double *priorgamscale, int varnum, gsl_matrix *posterior, double denominator,const int maxiters, const double epsabs);

int laplace_g_marg (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg_marg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg_marg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf_marg (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

void calc_Gaussiannode_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                         datamatrix *designmatrix, const double *priormean, const double *priorsd, const double *priorgamshape,
                          const double *priorgamscale, int varnum, int gaussiannodeid, gsl_matrix *posterior, double denominator,const int maxiters, const double epsabs);
                          
 int laplace_gaus_g_marg (const gsl_vector *beta, void *params,double *gvalue);

int laplace_gaus_dg_marg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_gaus_hessg_marg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_gaus_fdf_marg (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_gaus_inits_marg(gsl_vector *myBeta,struct fnparams *gparams); 
                          
