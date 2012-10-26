void calc_binary_marginal_rv_R(network *dag, datamatrix *obsdata, int nodeid,  int errverbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,double epsabs_inner, int maxiters_inner, double finitestepsize, int verbose,
				 double h_guess, double h_epsabs, int maxiters_hessian,
			       double *denom_modes, int paramid, double betafixed, double mlik, double *posterior);
				
double g_outer_marg_R(int Rn, double *betaincTauDBL, void *params);


void rv_dg_outer_marg_R(int n, double *betaDBL, double *dgvaluesDBL,void *params);

int rv_hessg_outer_marg( gsl_vector* betashort, void* params, gsl_matrix* hessgvalueshort,double h, gsl_matrix* hessgvalueshort3pt);

double compute_mlik_marg(double finitestepsize, void *params);

double compute_mlik_marg_mn(const gsl_vector *finitestepsize_vec, void *params);
