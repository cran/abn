#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
     
void calc_network_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix,  const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale,const int maxiters, const double epsabs, const int errverbose);


void calc_network_Score_laplace_reuse(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix,  const double *priormean, const double *priorsd, const double *priorgamshape, const double *priorgamscale, 
				 struct database *prevNodes,const int maxiters, const double epsabs, const int errverbose, int enforce_db_size);

double calc_node_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix,  const double *priormean, const double *priorsd,
                                 const double *priorgamshape, const double *priorgamscale,const int maxiters, const double epsabs);

void print_state (int iter, gsl_multiroot_fdfsolver * s);

int laplace_g (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

/*int generate_inits(gsl_vector *myBeta,datamatrix *designmatrix);*/

int generate_inits_n(gsl_vector *myBeta,struct fnparams *gparams);

double calc_Gaussiannode_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix,  const double *priormean, const double *priorsd,
                                 const double *priorgamshape, const double *priorgamscale, int gaussiannodeid,const int maxiters, const double epsabs);

 int laplace_gaus_g (const gsl_vector *beta, void *params,double *gvalue);

int laplace_gaus_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_gaus_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_gaus_fdf (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_gaus_inits(gsl_vector *myBeta,struct fnparams *gparams); 
