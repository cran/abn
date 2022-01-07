###############################################################################
## fitabn.R ---
## Author : Fraser Lewis & Gilles Kratzer
## Last modified : 02/12/2012
## Last modified : 29/09/2014 by Marta Pittavino, just renamed.
##               :  06/12/2016 by GK (implementation of formula
############################################################################### statment)

## fit a given DAG to data

fitabn <- function(...) {
    .Deprecated("fitAbn", msg="'fitabn' is deprecated.\n Use 'fitAbn' instead but note that arguments have slightly changed.")
    fitAbn(...)
}



fitAbn <- function(object=NULL, dag=NULL, data.df=NULL, data.dists=NULL, method=NULL, group.var=NULL, adj.vars=NULL, cor.vars=NULL,  centre=TRUE, create.graph=FALSE, compute.fixed=FALSE,
    control=list(), verbose=FALSE, ...) {


    ## method abnCache
    if (!is.null(object) & inherits(x=object,what="abnLearned")) {

        dag <- object$dag
        object <- object$score.cache

        data.df <- object$data.df
        data.dists <- object$data.dists
        # group.var <- object$group.var
        # cor.vars <- object$cor.vars
        # adj.vars <- object$adj.vars

        fitmethod <- object$method
    } else  fitmethod <- NULL

    if (is.null(method))  {
        method <- if (is.null(fitmethod)) "bayes" else fitmethod
    }else{
        if (verbose & (!is.null(fitmethod)))
            if (method != fitmethod) cat("Fitting and learned methods differ!\n")
    }


    ## start tests
    method <- c("bayes", "mle")[pmatch(tolower(method), c("bayes", "mle"))][1]
    if (is.na(method)) stop("'method' should be 'bayes' or 'mle'.")

    ## Control input
    if (method == "bayes") {
        con <- list(mean=0, prec=0.001, loggam.shape=1, loggam.inv.scale=5e-05, max.mode.error=10, max.iters=100, epsabs=1e-07, error.verbose=FALSE, epsabs.inner=1e-06, max.iters.inner=100,
            finite.step.size=1e-07, hessian.params=c(1e-04, 0.01), max.iters.hessian=10, max.hessian.error=1e-04, factor.brent=100, maxiters.hessian.brent=10, num.intervals.brent=100,
            min.pdf=0.001, n.grid=100, std.area=TRUE, marginal.quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), max.grid.iter=1000, marginal.node=NULL, marginal.param=NULL, variate.vec=NULL, seed=9062019)
    }

    if (method == "mle") {
        con <- list(maxit=100, tol=10^-11, seed=9062019)
    }

    nmsC <- names(con)

    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse=", "))


    if (method == "bayes") {

        out <- fitAbn.bayes(dag, data.df=data.df, data.dists=data.dists, group.var=group.var,
            cor.vars=cor.vars, create.graph=create.graph, compute.fixed=compute.fixed, mean=con$mean,
            prec=con$prec, loggam.shape=con$loggam.shape, loggam.inv.scale=con$loggam.inv.scale,
            max.mode.error=con$max.mode.error, max.iters=con$max.iters, epsabs=con$epsabs,
            error.verbose=con$error.verbose, epsabs.inner=con$epsabs.inner, max.iters.inner=con$max.iters.inner,
            finite.step.size=con$finite.step.size, hessian.params=con$hessian.params, max.iters.hessian=con$max.iters.hessian,
            max.hessian.error=con$max.hessian.error, factor.brent=con$factor.brent, maxiters.hessian.brent=con$maxiters.hessian.brent, num.intervals.brent=con$num.intervals.brent, min.pdf=con$min.pdf,
            n.grid=con$n.grid, std.area=con$std.area, marginal.quantiles=con$marginal.quantiles, max.grid.iter=con$max.grid.iter, marginal.node=con$marginal.node, marginal.param=con$marginal.param,
            variate.vec=con$variate.vec, seed=con$seed, verbose=verbose)

        return(out)
    }
    if (method == "mle") {

        out <- fitAbn.mle(dag, data.df=data.df, data.dists=data.dists, adj.vars=adj.vars,
                cor.vars=cor.vars, verbose=verbose, centre=centre, maxit=con$maxit, tol=con$tol, seed=con$seed)

        return(out)
    }

}  #EOF
