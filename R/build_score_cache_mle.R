############################################################################### buildscorecache.mle.R --- Author : Gilles Kratzer Last modified : 08/03/2017 Last modified : 18/03/2017 (mycache) Last modified : 20/10/2017 Multinomial

## fit a given DAG to data
buildscorecache.mle <- function(data.df = NULL, data.dists = NULL, max.parents = NULL, adj.vars = NULL, cor.vars = NULL, dag.banned = NULL, dag.retained = NULL, which.nodes = NULL, maxit = 100, 
    tol = 10^-8, centre = TRUE, defn.res = NULL, dry.run = FALSE, verbose = FALSE) {
    
    ## which.nodes
    if (!is.null(which.nodes)) {
        data.df <- data.df[, which.nodes]
        data.dists <- data.dists[which.nodes]
    }
    
    ## number of variable:
    n <- length(data.dists)
    nobs <- dim(data.df)[1]
    
    ## for compatibility
    group.var <- NULL
    
    # test same order for data.frame and data.dist
    if (Reduce("|", names(data.dists) != names(data.dists[colnames(data.df)]))) {
        stop("data.dists, data.df do not have the same names or the same names' order")
    }
    
    # formating factor
    data.df.lvl <- data.df
    
    ## standardize gaussian variables to zero mean and sd=1 have at least one gaussian variable
    if (centre && !is.null(data.dists == "gaussian")) {
        for (i in names(data.dists)[(data.dists == "gaussian")]) {
            data.df[, i] <- (data.df[, i] - mean(data.df[, i]))/sd(data.df[, i])
        }
    }
    
    for (i in 1:n) {
        if (data.dists[[i]] == "binomial" & class(data.df[, i]) != "numeric") {
            data.df[, i] <- as.numeric(factor(data.df[, i])) - 1
        }
        if (data.dists[[i]] == "multinomial") {
            data.df[, i] <- factor(data.df[, i])
        }
    }
    
    # adjustment: storing of df
    if (!is.null(adj.vars)) {
        data.df.adj <- data.df
        data.df <- data.df[, -adj.vars]
        n <- n - length(adj.vars)
    }
    
    ############## BAN / RETAIN
    
    # test for dag
    if (!is.null(dag.banned)) {
        if (is.matrix(dag.banned)) {
            ## run a series of checks on the DAG passed
            dag.banned <- check.valid.dag(dag.m = dag.banned, data.df = data.df, is.ban.matrix = TRUE, group.var = group.var)
        } else {
            if (grepl("~", as.character(dag.banned)[1], fixed = TRUE)) {
                dag.banned <- formula.abn(f = dag.banned, name = colnames(data.df))
                ## run a series of checks on the DAG passed
                dag.banned <- check.valid.dag(dag.m = dag.banned, data.df = data.df, is.ban.matrix = TRUE, group.var = group.var)
            }
        }
    } else {
        dag.banned <- check.valid.dag(dag.m = dag.banned, data.df = data.df, is.ban.matrix = TRUE, group.var = group.var)
    }
    
    # test for dag
    if (!is.null(dag.retained)) {
        if (is.matrix(dag.retained)) {
            ## run a series of checks on the DAG passed
            dag.retained <- check.valid.dag(dag.m = dag.retained, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var)
        } else {
            if (grepl("~", as.character(dag.retained)[1], fixed = T)) {
                dag.retained <- formula.abn(f = dag.retained, name = colnames(data.df))
                ## run a series of checks on the DAG passed
                dag.retained <- check.valid.dag(dag.m = dag.retained, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var)
            }
        }
    } else {
        dag.retained <- check.valid.dag(dag.m = dag.retained, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var)
    }
    
    ############################## Function to create the cache
    
    
    if (!is.null(defn.res)) {
        max.parents <- max(apply(defn.res[["node.defn"]], 1, sum))
        
    } else {
        
        ## max parents
        if (is.null(max.parents)) {
            max.parents <- 1
        }
        
        if (is.numeric(max.parents)) {
            if (max.parents >= n) {
                max.parents <- n - 1
            }
        }
        
        max.parent.list <- NULL
        
        if (is.list(max.parents)) {
            if (do.call(max, max.parents) > (n)) {
                stop("max.parent should be an integer or a list with a maximum possible of number of node-1 and the length of the list should not exceed the number of nodes.")
            } else {
                max.parent.list <- max.parents
                max.parents <- do.call(max, max.parents)
            }
        }
        
        
        ## Computing the cache
        fun.return <- function(x) {
            v <- rep(0, n - 1)
            v[x] <- 1
            return(v)
        }
        
        node.defn <- matrix(data = as.integer(0), nrow = 1L, ncol = n)
        children <- 1
        
        for (j in 1:n) {
            if (j != 1) {
                node.defn <- rbind(node.defn, matrix(data = as.integer(0), nrow = 1L, ncol = n))
                children <- cbind(children, j)
            }
            # node.defn <- rbind(node.defn,matrix(data = 0,nrow = 1,ncol = n))
            
            for (i in 1:(max.parents)) {
                
                tmp <- t(combn(n - 1, i, FUN = fun.return, simplify = TRUE))
                tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))
                
                node.defn <- rbind(node.defn, tmp)
                
                # children position
                children <- cbind(children, t(rep(j, length(tmp[, 1]))))
            }
        }
        # children <- rowSums(node.defn)
        colnames(node.defn) <- colnames(data.df)
        ## Coerce numeric matrix into integer matrix !!!
        node.defn <- apply(node.defn, c(1, 2), function(x) {
            (as.integer(x))
        })
        
        children <- as.integer(children)
        # node.defn_ <- node.defn
        
        ## DAG RETAIN/BANNED
        for (i in 1:n) {
            for (j in 1:n) {
                
                ## DAG RETAIN
                if (dag.retained[i, j] != 0) {
                  tmp.indices <- which(children == i & node.defn[, j] == 0)
                  
                  if (length(tmp.indices) != 0) {
                    node.defn <- node.defn[-tmp.indices, ]
                    children <- children[-tmp.indices]
                  }
                }
                
                ## DAG BANNED
                if (dag.banned[i, j] != 0) {
                  tmp.indices <- which(children == i & node.defn[, j] == 1)
                  
                  if (length(tmp.indices) != 0) {
                    node.defn <- node.defn[-tmp.indices, ]
                    children <- children[-tmp.indices]
                  }
                }
                
            }
        }
        
        mycache <- list(children = as.integer(children), node.defn = (node.defn))
        
        ###------------------------------###
        ### start limiting max.parent list###
        ###------------------------------###
        
        if (!is.null(max.parent.list)) {
            for (z in 1:n) {
                tmp <- mycache[["node.defn"]][mycache[["children"]] == z, ]
                
                mycache[["node.defn"]][mycache[["children"]] == z, ] <- tmp[rowSums(tmp) <= unlist(max.parent.list[z]), ]
            }
        }
        
        
        ###----------------###
        ### start adjustment###
        ###----------------###
        
        if (!is.null(adj.vars)) {
            
            # mycache$node.defn.adj <- mycache$node.defn
            
            ## adding adjustment column set to zero mycache$node.defn.adj <- cbind(mycache$node.defn,matrix(data = 0,nrow = dim(mycache$node.defn)[1],ncol = length(adj.vars)))
            mycache$node.defn <- cbind(mycache$node.defn, matrix(data = 0, nrow = dim(mycache$node.defn)[1], ncol = length(adj.vars)))
            
            if (is.null(cor.vars)) {
                cor.vars <- colnames(data.df)
            }
            
            ## adjustment variables
            
            mycache$node.defn[mycache$children[match(cor.vars, colnames(data.df))], dim(data.df)[2]:dim(data.df.adj)] <- 1
            
            ## output
            colnames(mycache$node.defn) <- c(colnames(data.df), adj.vars)
            
            mycache$node.defn <- mycache$node.defn[, names(data.df.adj)]
            data.df <- data.df.adj
        }
        
        ##----------------------
        ## multinomial adaptation
        ##----------------------
        
        # unpacking the multinomial variables in the cache
        repetition.multi <- vector(length = n)
        
        for (i in 1:n) {
            if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
                repetition.multi[i] <- 1
            } else {
                repetition.multi[i] <- nlevels(data.df.lvl[, i])
            }
        }
        
        if (!is.null(adj.vars)) {
            mycache$node.defn.multi <- mycache$node.defn.adj[, rep(1:n, repetition.multi)]
            data.df <- data.df.adj[, colnames(mycache$node.defn.adj)]
        } else {
            mycache$node.defn.multi <- mycache$node.defn[, rep(1:n, repetition.multi)]
            
        }
        
        # unpacking the multinomial variables in the data.df
        
        data.df.multi <- NULL
        
        for (i in 1:n) {
            if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
                data.df.multi <- as.data.frame(cbind(data.df.multi, data.df[, i]))
                colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
            } else {
                tmp <- model.matrix(~-1 + factor(data.df.lvl[, i]))
                colnames(tmp) <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
                data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))
            }
        }
        
    }
    if (dry.run) {
        return(mycache)
    }
    
    ## EOF cache creation
    
    out <- list()
    ll.out.tmp <- list()
    aic.out.tmp <- list()
    bic.out.tmp <- list()
    mdl.out.tmp <- list()
    row.num.out <- list()
    child.out <- list()
    
    row.num <- 1
    
    for (child in mycache[["children"]]) {
        
        Y <- data.matrix(data.df[, child])
        
        if (is.null(adj.vars)) {
            if ("multinomial" %in% data.dists[as.logical(mycache$node.defn[row.num, ])]) {
                X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
            } else {
                X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
            }
        }
        
        if (!is.null(adj.vars)) {
            if ("multinomial" %in% data.dists[as.logical(mycache$node.defn.adj[row.num, ])]) {
                X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
            } else {
                X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
            }
        }
        
        ## Rank deficiency
        
        num.na <- 0
        
        if (qr(X)$rank/ncol(X) != 1 & as.character(data.dists[child]) == "binomial") {
            
            Y <- as.numeric(as.character(Y))
            
            tryCatch(fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = maxit, tol = tol), error = function(e) {
                while ((qr(X)$rank/ncol(X)) != 1) {
                  X <- X[, -1]
                  num.na <- num.na + 1
                  if (is.null(ncol(X))) 
                    X <- as.matrix(X)
                }
                fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = maxit, tol = tol)
            }, finally = fit)
            
        } else {
            
            switch(as.character(data.dists[child]), binomial = {
                
                Y <- as.numeric(as.character(Y))
                
                fit <- tryCatch({
                  (fit <- irls_binomial_cpp_fast(A = X, b = Y, maxit = maxit, tol = tol))
                }, error = function(e) {
                  fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = maxit, tol = tol)
                  return(fit)
                })
                
                if (is.na(sum(fit[[1]]))) fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = maxit, tol = tol)
                
            }, gaussian = {
                
                suppressWarnings(fit <- irls_gaussian_cpp_fast(A = X, b = Y, maxit = maxit, tol = tol))
                
            }, poisson = {
                
                suppressWarnings(fit <- irls_poisson_cpp_fast(A = X, b = Y, maxit = maxit, tol = tol))
                
            }, multinomial = {
                
                Ymulti <- data.matrix(model.matrix(~-1 + data.df.lvl[, child]))
                
                # weights w <- rep(1, dim(Ymulti)[1L])
                
                r <- ncol(X)
                p <- ncol(Ymulti)
                mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE, r)), p - 1L))
                
                tmp <- nnet.default(x = X, y = Ymulti, mask = mask, size = 0, skip = TRUE, softmax = TRUE, rang = 0, trace = FALSE)
                
                fit <- NULL
                
                fit$loglik <- -tmp$value
                
                edf <- ifelse(length(tmp$lev) == 2L, 1, length(tmp$lev) - 1) * qr(X)$rank
                
                fit$aic <- 2 * tmp$value + 2 * edf
                
                fit$bic <- 2 * tmp$value + edf * log(nobs)
                
            })
        }
        
        ll.out.tmp[[row.num]] <- fit$loglik
        aic.out.tmp[[row.num]] <- fit$aic
        bic.out.tmp[[row.num]] <- fit$bic
        mdl.out.tmp[[row.num]] <- fit$bic + (1 + sum(mycache[["node.defn.multi"]][row.num, ]) - num.na) * log(n)
        row.num.out[[row.num]] <- row.num
        child.out[[row.num]] <- child
        
        row.num <- row.num + 1
        
        ## reinitialisation
        num.na <- 0
    }
    
    out[["children"]] <- as.integer(child.out)
    out[["node.defn"]] <- mycache$node.defn
    out[["mlik"]] <- as.numeric(ll.out.tmp)
    out[["error.code"]] <- list()
    out[["hessian.accuracy"]] <- list()
    out[["used.INLA"]] <- list()
    out[["error.code.desc"]] <- list()
    out[["data.df"]] <- data.df.lvl
    out[["data.dists"]] <- data.dists
    out[["max.parents"]] <- max.parents
    out[["dag.retained"]] <- dag.retained
    out[["dag.banned"]] <- dag.banned
    out[["group.ids"]] <- list()
    out[["aic"]] <- as.numeric(aic.out.tmp)
    out[["bic"]] <- as.numeric(bic.out.tmp)
    out[["mdl"]] <- as.numeric(mdl.out.tmp)
    
    out[["method"]] <- "mle"
    
    return(out)
    
}

