void init_nodedatabase(struct database *prevNodes,network *dag, int  db_size, int allocMem);
void storenodescore(network *dag, int nodeid,double nodescore, struct database *prevNodes, int enforce_db_size);
int nodescoreisknown(network *dag, int nodeid,double *nodescore, struct database *prevNodes);
