############################################################################### markov-blanket.R --- Author : Gilles Kratzer Document created : 12/04/2017 Last modification : 12/04/2017


##---------------------------------------------------------------------------
## Function that return the markov blanket of a given node in a given network
##---------------------------------------------------------------------------


mb <- function(dag, node, data.dists=NULL) {
    
    dag <- validate_abnDag(dag,  data.df=data.dists, returnDAG=TRUE)

    if (!all(node%in%colnames(dag)))
        stop("Incorrect node specification 'node' to compute the Markov Blanket")
    
    ## as of older version.
    # row parent column children mb.node.final <- list()
    mb.node.tmp <- list()
    for (n.element in node) {
        
        ## Parent + Children
        mb.children <- list()
        mb.parent <- list()
        for (i in 1:length(dag[1, ])) {
            if (dag[i, n.element] != 0) {
                mb.children[i] <- names(dag[, n.element])[i]
            }
            if (dag[n.element, i] != 0) {
                mb.parent[i] <- names(dag[n.element, ])[i]
            }
        }
        # delete NULL element
        mb.children <- unlist(mb.children[!sapply(mb.children, is.null)])
        mb.parent <- unlist(mb.parent[!sapply(mb.parent, is.null)])
        
        ## Parent of children
        mb.parent.children <- list()
        for (node.children in mb.children) {
            for (i in 1:length(dag[1, ])) {
                if (dag[node.children, i] != 0) {
                    mb.parent.children[i] <- names(dag[node.children, ])[i]
                }
            }
        }
        # delete NULL element
        mb.parent.children <- unlist(mb.parent.children[!sapply(mb.parent.children, is.null)])
        
        # add all list
        mb.node <- unlist(list(mb.children, mb.parent, mb.parent.children))
        
        # unique element
        mb.node <- unique(mb.node)
        
        # delete index node
        mb.node.wo <- NULL
        if (length(mb.node) != 0) {
            for (i in 1:length(mb.node)) {
                if (mb.node[c(i)] == n.element) {
                    mb.node.wo <- mb.node[-c(i)]
                }
            }
        }
        if (is.null(mb.node.wo)) {
            mb.node.wo <- mb.node
        }
        
        ## store out of loop
        mb.node.tmp <- unlist(list(mb.node.tmp, mb.node.wo))
        # unique element
        mb.node.tmp <- unique(mb.node.tmp)
        # delete NULL element
        mb.node.tmp <- unlist(mb.node.tmp[!sapply(mb.node.tmp, is.null)])
        
    }  #EOF loop through node
    
    return(mb.node.tmp)
}


# EOF
