
void hill_climb_iter_additive(storage *nodescore,cycle *cyclestore,network *dag_orig, network *dag_scratch,network *dag_opt1,network *dag_opt2, network *dag_opt3, unsigned                               int maxparents, datamatrix *obsdata, datamatrix *designmatrix, const double *priormean, const double *priorsd,
                              const double *priorgamshape, const double *priorgamscale, int verbose, struct database *prevNodes,
                              const int maxiters, const double epsabs, const int errverbose, int enforce_db_size);
