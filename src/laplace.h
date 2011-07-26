#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
     
void calc_network_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd);

double calc_node_Score_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd);

void print_state (int iter, gsl_multiroot_fdfsolver * s);

int laplace_g (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_inits(gsl_vector *myBeta, datamatrix *designmatrix);
