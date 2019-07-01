############################################################################### 
## abn-toolbox.R --- 
## Author : Gilles Kratzer 
## Document created : 01/02/2017 
## Last modification : 01/02/2017 
## Last modification : 11/04/2017 (compareDag, infoDag) 
## Last modification : 21/04/2017 (simulateDag) 
## Last modification : 29/08/2017 (compareDag (G and F1 score), essential.graph (minimal vs completed)) 
## Last modification : 16/07/2018 error in Hamming distance corrected
############################################################################### 

##-------------------------------------------------------------------------
## External usefull functions to analayse ABN and BN
##-------------------------------------------------------------------------

# functions logit
logit <- function(x) {
    log(x/(1 - x))
}

# functions expit
expit <- function(x) {
    exp(x)/(1 + exp(x))
}

# function that compute an odds ratio from a table/matrix
or <- function(x) {
    
    x <- as.matrix(x)
    
    if (dim(x)[1] != 2 || dim(x)[2] != 2) {
        stop("The contengency table should be of 2-dimensionnal with a 2x2 formulation.")
    }
    
    if (x[1, 2] == 0) {
        x[1, 2] <- 0.5
    }
    if (x[2, 1] == 0) {
        x[2, 1] <- 0.5
    }
    
    out <- (x[1, 1] * x[2, 2])/(x[1, 2] * x[2, 1])
    
    return(out)
}

# Probability to odds
odds <- function(x) {
    (x/(1 - x))
}

##-------------------------------------------------------------------------
## External function that compares DAGs and computes many metrics A Comparison of Structural Distance Measures for Causal Bayesian Network Models
##-------------------------------------------------------------------------


compareDag <- function(ref, test, node.names = NULL) {
    
    ## check ref dag
    ref <- validate_abnDag(  ref, data.df=node.names, returnDAG=TRUE)
    test <- validate_abnDag( test, data.df=node.names, returnDAG=TRUE)

    if (any(dim(ref) != dim(test))) {
        stop("The reference or test DAG has not the same size")
    }
    
    n <- dim(ref)[1]
    
    ## unit matrix
    ref[ref != 0] <- 1
    test[test != 0] <- 1
    
    diff.matrix <- ref - (0.5 * test)
    
    diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))
    
    if(sum(ref == 1)==0 | sum(test == 1)==0){
        warning("If the test or reference matrix is an empty matrix some of the estimates are not defined.")
    }
    
    ## output
    out <- list()
    
    out[["TPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1))
    out[["FPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == -0.5]))/(sum(ref == 0))
    out[["Accuracy"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]) + as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0]))/(dim(ref)[1]^2)
    out[["FDR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))
    out[["G-measure"]] <- sqrt(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1)) * (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))
    out[["F1-score"]] <- (2/((1/as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))) + (1/(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))))
    out[["PPV"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))
    out[["FOR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))
    
    # transforming "reverse arc" into -0.5
    for( i in 1:n){
        for( j in i:n){
            if( diff.matrix[i,j]!=0 & diff.matrix[j,i]!=0){
                diff.matrix[i,j] <- -(diff.matrix[i,j] + diff.matrix[j,i])
                diff.matrix[j,i] <- 0
            }
        }
    }
    
    diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))
    
    out[["Hamming-distance"]] <- sum(as.numeric(diff.matrix.tab[names(diff.matrix.tab) %in% c(-0.5, 1)]))
    
    return(out)
}


##-------------------------------------------------------------------------
## External function that return information about a dag number of nodes number of arcs avergae size of the MB average size of the Neighborhood
##-------------------------------------------------------------------------

infoDag <- function(dag, node.names = NULL) {
    
    ## dag transformation
    if (!is.null(dag)) {
        if (is.matrix(dag)) {
            ## run a series of checks on the DAG passed
            dag <- abs(dag)
            diag(dag) <- 0
            dag <- check.valid.dag(dag.m = dag, is.ban.matrix = FALSE, group.var = NULL)
            ## naming
            if (is.null(colnames(dag))) {
                colnames(dag) <- rownames(dag) <- node.names
            }
        } else {
            if (grepl("~", as.character(dag)[1], fixed = TRUE)) {
                dag <- formula.abn(f = dag, name = node.names)
                ## run a series of checks on the DAG passed
                dag <- check.valid.dag(dag.m = dag, is.ban.matrix = FALSE, group.var = NULL)
            }
        }
    } else {
        stop("Dag specification must either be a matrix or a formula expression")
    }
    
    dag[dag != 0] <- 1
    diag(dag) <- 0
    
    if (is.null(node.names) & is.null(colnames(dag))) {
        stop("Either name have to be provided with a formula statement either a named matrix to define a DAG")
    }
    if (is.null(node.names)) {
        node.names <- colnames(dag)
    }
    
    ## ======================== test for conformity over! ========================
    
    out <- list()
    ## number of nodes
    n.nodes <- dim(dag)[1]
    
    ## number of arcs
    n.arcs <- sum(dag)
    
    ## =========================== average markov blanket size ===========================
    
    mb.size <- vector(mode = "numeric", length = length(node.names))
    
    for (i in 1:length(node.names)) {
        
        node <- node.names[i]
        # row children column parent
        
        ## Parent + Children
        mb.children <- list()
        mb.parent <- list()
        for (j in 1:length(dag[1, ])) {
            if (dag[node, j] != 0) {
                mb.children[j] <- names(dag[node, ])[j]
            }
            if (dag[j, node] != 0) {
                mb.parent[j] <- names(dag[, node])[j]
            }
        }
        # delete NULL element
        mb.children <- unlist(mb.children[!sapply(mb.children, is.null)])
        mb.parent <- unlist(mb.parent[!sapply(mb.parent, is.null)])
        
        ## Parent of children
        mb.parent.children <- list()
        for (node.children in mb.children) {
            for (k in 1:length(dag[1, ])) {
                if (dag[k, node.children] != 0) {
                  mb.parent.children[k] <- names(dag[, node.children])[k]
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
            for (l in 1:length(mb.node)) {
                if (mb.node[c(l)] == node) {
                  mb.node.wo <- mb.node[-c(l)]
                }
            }
        }
        if (is.null(mb.node.wo)) {
            mb.node.wo <- mb.node
        }
        
        mb.size[i] <- length(mb.node.wo)
    }
    
    mb.average <- mean(mb.size)
    
    ## average Neighborhood
    
    nh.size <- vector(mode = "numeric", length = length(node.names))
    parent.size <- vector(mode = "numeric", length = length(node.names))
    children.size <- vector(mode = "numeric", length = length(node.names))
    
    for (i in 1:length(node.names)) {
        nh.size[i] <- sum(dag[i, ]) + sum(dag[, i])
        parent.size[i] <- sum(dag[i, ])
        children.size[i] <- sum(dag[, i])
    }
    nh.average <- mean(nh.size)
    parent.average <- mean(parent.size)
    children.average <- mean(children.size)
    
    ## output
    out[["n.nodes"]] <- n.nodes
    out[["n.arcs"]] <- n.arcs
    out[["mb.average"]] <- mb.average
    out[["nh.average"]] <- nh.average
    out[["parent.average"]] <- parent.average
    out[["children.average"]] <- children.average
    
    return(out)
}

##-------------------------------------------------------------------------
## External function that simulate a dag with arbitrary arcs density
##-------------------------------------------------------------------------

simulateDag <- function(node.name = NULL, data.dists = NULL, nn = 0.5) {
    
    ## test
    if (length(node.name) <= 2) {
        stop("No need for help to simulate a DAG of size 2")
    }
    if (nn > 1 | nn < 0) {
        stop("The network density should be a real number in [0,1]")
    }
    
    if (is.null(data.dists)) {
        data.dists <- sample(x = c("gaussian", "binomial", "poisson"), size = length(node.name), replace = TRUE)
        names(data.dists) <- node.name
    }
    
    dag <- matrix(data = 0, nrow = length(node.name), ncol = length(node.name))
    
    
    for (j in 2:(length(node.name))) {
        dag[c(j:length(node.name)), (j - 1)] <- rbinom(n = (length(node.name) - j + 1), size = 1, prob = nn)
    }
    order.m <- sample(x = length(node.name), size = length(node.name), replace = FALSE)
    
    ## changing order
    dag <- dag[order.m, order.m]  #order.m
    
    ## naming
    colnames(dag) <- rownames(dag) <- node.name
    
    out <- create_abnDag(dag, data.dists=data.dists)
    ## structure
    return(out)
}

##-------------------------------------------------------------------------
## Function that computes skewness of a distribution
##-------------------------------------------------------------------------

skewness <- function(x) {
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
}

##-------------------------------------------------------------------------
## External function that computes essential graph of a dag Minimal PDAG: The only directed edges are those who participate in v-structure Completed PDAG: very directed edge corresponds to a
## compelled edge, and every undirected edge corresponds to a reversible edge
##-------------------------------------------------------------------------

essentialGraph <- function(dag, node.names = NULL, PDAG = "minimal") {

    ## dag transformation
    if (is.matrix(dag)) {
        # check.valid.dag(dag.m=dag.m,data.df=data.df,is.ban.matrix=FALSE,group.var=group.var) unit matrix
        dag[dag != 0] <- 1
        node.names <- colnames(dag)
        
    } else {
        if (grepl("~", as.character(dag)[1], fixed = TRUE)) {
            dag <- formula.abn(f = dag, name = node.names)
        } else {
            stop("Dag specification must either be a matrix or a formula expresion")
        }
    }
    
    ## compute essential graph moral graph
    
    dim.dag <- dim(dag)[1]
    moral <- matrix(data = 0, nrow = dim.dag, ncol = dim.dag)
    
    for (i in 1:dim.dag) {
        for (j in 1:dim.dag) {
            if (dag[i, j] == 1) {
                moral[i, j] <- 1
                moral[j, i] <- 1
            }
        }
    }
    
    colnames(moral) <- rownames(moral) <- node.names
    
    ## essential arcs
    if (PDAG == "minimal") {
        for (i in 1:dim.dag) {
            if (sum(dag[i, ]) >= 2) {
                for (j in 1:dim.dag) {
                  if (dag[i, j] == 1) {
                    moral[j, i] <- 0
                  }
                }
            }
        }
        
        
        colnames(moral) <- rownames(moral) <- node.names
        return(moral)
    }
    
    if (PDAG == "completed") {
        for (i in 1:dim.dag) {
            if (sum(dag[i, ]) >= 2) {
                for (j in 1:dim.dag) {
                  if (dag[i, j] == 1) {
                    moral[j, i] <- 0
                  }
                  if (dag[j, i] == 1) {
                    moral[i, j] <- 0
                  }
                }
            }
        }
        
        colnames(moral) <- rownames(moral) <- node.names
        return(moral)
    }
}

##-------------------------------------------------------------------------
##Function that compute likelihood contribution of observation
##-------------------------------------------------------------------------

scoreContribution <- function(object = NULL, 
                              dag.m = NULL, 
                              data.df = NULL,
                              data.dists = NULL,
                              verbose = FALSE){
    
    ## method abnCache
    if (!is.null(object)){
        if (inherits(object, "abnLearned")) {
            dag.m <- object$dag
            data.df <- object$score.cache$data.df
            data.dists <- object$score.cache$data.dists
        }}
    
    ## transform factor into 0/1
    
    node.ordering <- names(data.dists)
    nobs <- dim(data.df)[1]
    
    ll <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    nb.param <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    nb.parents <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    colnames(ll) <- colnames(nb.param) <- colnames(nb.parents) <- node.ordering
    nb.parents <- rowSums(dag.m)/nobs
    for (node in node.ordering) {
        
        if(as.character(data.dists[node])=="binomial"){
            Y <- as.factor(data.matrix(data.df[, node]))
        } else {
            Y <- data.matrix(data.df[, node])
        }
        X <- data.matrix(cbind(rep(1,nobs),data.df[, as.logical(dag.m[node,])]))
        if(as.character(data.dists[node])=="gaussian"){
            nb.param[,node] <- (dim(X)[2] + 1)/nobs
        }else{
            nb.param[,node] <- dim(X)[2]/nobs
        }
        
        x <- glm(formula = Y ~ -1 + X,family = as.character(data.dists[node]))
        yhat <- predict.glm(x,newdata = x$data, type = "response")
        
        switch (as.character(data.dists[node]),
                "binomial" = {
                    ll[,node] <- dbinom(as.numeric(Y)-1, size = 1L, prob = yhat, log = TRUE)
                },
                "gaussian" = {
                    ll[,node] <- dnorm(Y, mean = yhat, sd = sigma(x),log = TRUE)              
                },
                "poisson" = {
                    ll[,node] <- dpois(Y, lambda = yhat, log = TRUE)
                }
        )
        hv <- hatvalues(x)
        
    }
    
    aic <- - 2*ll + 2*nb.param
    bic <- - 2*ll + nb.param*(log(nobs))
    mdl <- bic + (1/nobs + nb.parents) * log(nobs)
    
    out <- list("mlik" = ll, "aic" = aic, "bic" = bic, "mdl" = mdl, "hatvalues" = hv)
    
    return(out)
    
}


##-------------------------------------------------------------------------
## Function to export abn coefficients to xls
##-------------------------------------------------------------------------



## EOF
