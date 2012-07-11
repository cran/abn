
double qprime(int *parents, int numNodes, double offset);

double feature(int *parents, int nodeid, int child, int parent, int i, int *customfeature, int numRows);

void mobius_transform(double *f_hat_0,double *f_hat_1, int **parentsloc, double *score, int numsubsets, int *nodeid, int numnodes,double **alpha);

int indexjcomplement(int currow, int curparent, int **parents, int numnodes, int numsubsets);

double g(int *V, int len, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, double *ptr_score);

#ifdef CHECKRECURSION
double c(int dropped, int *V, int *V2, int lenV);
#endif

double get_q(int dropped, int *V, int *V2, int lenV);

double get_alpha(int dropped, int *V, int *V2, int lenV, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, int *found);
#ifdef JUNK2
double get_alpha_maxparent_exceeded(int dropped, int *V, int *V2, int lenV, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents);
#endif
void mobius_transform_max(double *f_hat_0,double *f_hat_1, int **parentsloc, double *score, int numsubsets, int *nodeid, int numnodes, double **alpha);

double g_max(int *V, int len, double **alpha, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, int previousg, int toplevellen, int *maxnode, double *ptr_score);

int get_alpha_parents(int dropped, int *V, int *V2, int lenV, int **alphalookup, int **parents, int *parents_numparents, int *nodesid, int numRows, int numNodes, int maxparents, double *ptr_score);

double get_alpha_no_mobius(int node, int *V2, int **parents, int *parents_numparents, int *nodesid, double *ptr_score, int numRows,int numNodes,int maxparents, int lenV2);

double get_alpha_no_mobius_max(int node, int *V2, int **parents, int *parents_numparents, int *nodesid, double *ptr_score, int numRows,int numNodes,int maxparents, int lenV2);


