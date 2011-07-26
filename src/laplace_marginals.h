#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
     
void calc_network_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata, int verbose,
				 datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd,
				 int nodenum, int varnum, gsl_matrix *posterior);

void calc_node_Marginals_laplace(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid, int verbose,
                                datamatrix *designmatrix, SEXP R_labels, const double *priormean, const double *priorsd, int varnum, gsl_matrix *posterior, double denominator);

int laplace_g_marg (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg_marg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg_marg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf_marg (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

