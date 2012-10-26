void calc_network_Score(network *dag, datamatrix *obsdata, datamatrix *designmatrix,
				const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs, int verbose, const int errverbose, SEXP results, int storeModes, 
			double epsabs_inner, int maxiters_inner, double finitestepsize,double h_guess, double h_epsabs, int max_iters_hessian);
				
void calc_parameter_marginal(network *dag,datamatrix *obsdata, datamatrix *designmatrix,
				const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs, int verbose, const int errverbose, 
			      double *denom_modes, int childid, int paramid,  
			     double epsabs_inner, int maxiters_inner, double finitestepsize,
			      double h_guess, double h_epsabs, int maxiters_hessian,
			     double betafixed, double mlik, double *posterior);				
void storeResults(SEXP results,network *dag,int storeModes, int nodeid, int vartype);				
				
				
	 

	 
