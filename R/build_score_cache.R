############################################################################### buildscorecache.R --- Author : Gilles Kratzer Last modified : 21/05/2019

## precompute a cache of scores to data
buildscorecache <- function(data.df = NULL, data.dists = NULL, method = "bayes", group.var = NULL, adj.vars = NULL, cor.vars = NULL, dag.banned = NULL, dag.retained = NULL, max.parents = NULL, which.nodes = NULL, 
    defn.res = NULL, centre = TRUE, dry.run = FALSE, control = list(), verbose = FALSE, ...) {
    
    ## start tests
    method <- tolower(method)
    method <- c("bayes", "mle")[pmatch(method, c("bayes", "mle"))]
    
    data.dists <- validate_dists( data.dists, returnDists=TRUE)
    
    if (method == "bayes") {
        con <- list(max.mode.error = 10, mean = 0, prec = 0.001, loggam.shape = 1, loggam.inv.scale = 5e-05, max.iters = 100, epsabs = 1e-07, error.verbose = FALSE, epsabs.inner = 1e-06, max.iters.inner = 100, 
            finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01), max.iters.hessian = 10, max.hessian.error = 0.5, factor.brent = 100, maxiters.hessian.brent = 100, num.intervals.brent = 100)
    }
    
    if (method == "mle") {
        con <- list(maxit = 100, tol = 10^-8)
    }
    
    nmsC <- names(con)
    
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    
    
    if (dry.run == TRUE) {
        # return: list of child and parent combinations
        
        if (method == "bayes") {
            
            out <- buildscorecache.bayes(data.df = data.df, data.dists = data.dists, group.var = group.var, cor.vars = cor.vars, dag.banned = dag.banned, dag.retained = dag.retained, max.parents = max.parents, 
                which.nodes = which.nodes, defn.res = defn.res, dry.run = dry.run, max.mode.error = con$max.mode.error, verbose = verbose, centre = centre, mean = con$mean, prec = con$prec, loggam.shape = con$loggam.shape, 
                loggam.inv.scale = con$loggam.inv.scale, max.iters = con$max.iters, epsabs = con$epsabs, error.verbose = con$error.verbose, epsabs.inner = con$epsabs.inner, max.iters.inner = con$max.iters.inner, 
                finite.step.size = con$finite.step.size, hessian.params = con$hessian.params, max.iters.hessian = con$max.iters.hessian, max.hessian.error = con$max.hessian.error, factor.brent = con$factor.brent, 
                maxiters.hessian.brent = con$maxiters.hessian.brent, num.intervals.brent = con$num.intervals.brent)
            return(out)
        }
        if (method == "mle") {
            
            out <- buildscorecache.mle(data.df = data.df, data.dists = data.dists, adj.vars = adj.vars, cor.vars = cor.vars, dag.banned = dag.banned, dag.retained = dag.retained, max.parents = max.parents, 
                which.nodes = which.nodes, defn.res = defn.res, dry.run = dry.run, verbose = verbose, centre = centre, maxit = con$maxit, tol = con$tol)
            return(out)
        }
        
    } else {
        # compute scores
        
        if (method == "bayes") {
            
            out <- buildscorecache.bayes(data.df = data.df, data.dists = data.dists, group.var = group.var, cor.vars = cor.vars, dag.banned = dag.banned, dag.retained = dag.retained, max.parents = max.parents, 
                which.nodes = which.nodes, defn.res = defn.res, dry.run = dry.run, max.mode.error = con$max.mode.error, verbose = verbose, centre = centre, mean = con$mean, prec = con$prec, loggam.shape = con$loggam.shape, 
                loggam.inv.scale = con$loggam.inv.scale, max.iters = con$max.iters, epsabs = con$epsabs, error.verbose = con$error.verbose, epsabs.inner = con$epsabs.inner, max.iters.inner = con$max.iters.inner, 
                finite.step.size = con$finite.step.size, hessian.params = con$hessian.params, max.iters.hessian = con$max.iters.hessian, max.hessian.error = con$max.hessian.error, factor.brent = con$factor.brent, 
                maxiters.hessian.brent = con$maxiters.hessian.brent, num.intervals.brent = con$num.intervals.brent)
            
            class(out) <- c("abnCache")
            return(out)
        }
        if (method == "mle") {
            
            out <- buildscorecache.mle(data.df = data.df, data.dists = data.dists, adj.vars = adj.vars, cor.vars = cor.vars, dag.banned = dag.banned, dag.retained = dag.retained, max.parents = max.parents, 
                which.nodes = which.nodes, defn.res = defn.res, dry.run = dry.run, verbose = verbose, centre = centre, maxit = con$maxit, tol = con$tol)
            
            class(out) <- c("abnCache")
            return(out)
        }
    }
    
}  #EOF
