
void hill_climb_iter_additive(storage *nodescore,cycle *cyclestore,network *dag_orig, network *dag_scratch,network *dag_opt1,network *dag_opt2, network *dag_opt3, unsigned int maxparents, 
                    datamatrix *obsdata, datamatrix *designmatrix, const double *priormean, const double *priorsd, int verbose, SEXP R_labels);
