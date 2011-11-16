
void df_to_dm(SEXP R_obsdata,datamatrix *obsdata, SEXP R_numVarLevels);
void store_results(SEXP R_listresults,network *dag_best, int iter, SEXP ans, int verbose);
void df_to_dm_mixed(SEXP R_obsdata,datamatrix *obsdata, const int *vartype);

