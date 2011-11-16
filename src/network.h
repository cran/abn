
void build_init_dag(network *dag,datamatrix *obsdata,int maxparents);
void get_dag(network *dag,SEXP R_dag);
void get_dag_list(network *dag,SEXP R_dag,unsigned int index);
unsigned int hascycle(cycle *cyclestore,network *dag);
void get_numincomingedges(unsigned int *incomingedges,unsigned int **graph, unsigned int numnodes);
void calc_network_Score(storage *nodescore, network *dag,datamatrix *obsdata, double priordatapernode, int useK2, 
                         int verbose, SEXP R_labels);
double calc_node_Score(storage *nodescore,network *dag,datamatrix *obsdata,int nodeid,double priordatapernode, int useK2, 
                         int verbose, SEXP R_labels);
int crossproduct(storage *nodescore,int *numlevels,int **dest, int n,int maxparents);
void setbanlist(network *dag,SEXP R_dag);
void generate_random_dag(cycle *cyclestore,storage *nodescore,network *dag,unsigned nopermuts, unsigned int maxparents, SEXP R_shuffle, unsigned int offset,
                         unsigned int searchnum, SEXP R_dag_startlist);
void nullnetworkdefn(network *dag,cycle *cyclestore, unsigned int searchnum, SEXP R_dag_startlist);
unsigned int add_arc(cycle *cyclestore,network *dag, unsigned int child, unsigned int parent, unsigned int maxparents);
void hill_climb_iter(storage *nodescore,cycle *cyclestore,network *dag, network *dag_scratch,network *dag_opt1,network *dag_op2, network *dag_opt3, 
		     unsigned int maxparents, datamatrix *obsdata, double priordatapernode, int useK2, int verbose, SEXP R_labels, struct database *prevNodes, int enforce_db_size);
unsigned int remove_arc(network *dag, unsigned int child, unsigned int parent);
unsigned int reverse_arc(cycle *cyclestore,network *dag, unsigned int child, unsigned int parent, unsigned int maxparents);

void copynetworkdefn(network *dag,network *dag_scratch);

void init_hascycle(cycle *cyclestore, network *dag);
void init_network_score(storage *nodescore,network *dag);
void init_random_dag(storage *nodescore,network *dag);
void setretainlist(network *dag,SEXP R_dag);
void setstartlist(network *dag,SEXP R_dag);
void calc_network_Score_reuse(storage *nodescore, network *dag,datamatrix *obsdata, double priordatapernode, int useK2, 
                         int verbose, SEXP R_labels, struct database *prevNodes, int enforce_db_size);


