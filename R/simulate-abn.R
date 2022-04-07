################################################################################ simulate-abn.R --- Author : Gilles Kratzer Last Modified on: 06/12/2016

##-------------------------------------------------------------------------
## Function that simulate a dag from a matrix or a formula
##-------------------------------------------------------------------------

simulateAbn <- function(data.dists = NULL, data.param = NULL, data.param.var = NULL,
                        data.param.mult = NULL, n.chains = 10, n.adapt = 1000, n.thin = 100,
                        n.iter = 10000, bug.file = NULL, verbose = TRUE, simulate = TRUE, keep.file=FALSE, seed = 42) {

    if (!requireNamespace("rjags", quietly = TRUE)) {
        stop("library rjags is not available!\n")
    }

    ##-----------------------------------------------------------------------------------------
    ## initialisation
    ##-----------------------------------------------------------------------------------------

    ## for testing
    group.var <- NULL

    ## creating dag from data.param

    dag.m <- as.matrix(data.param)
    if (any(is.na(dag.m))) {
        stop("Data param specification must be a matrix that contain no Na's")
    }
    dag.m[dag.m != 0] <- 1
    if (is.null(colnames(dag.m))) {
        colnames(dag.m) <- names(data.dists)
        rownames(dag.m) <- names(data.dists)
    }
    diag(dag.m) <- 0

    ## check dag specifications
    dag.m <- check.valid.dag(dag = dag.m, is.ban.matrix = FALSE, group.var = group.var)

    # Missing precision parameter
    if (missing(data.param.var)) {
        data.param.var <- matrix(data = 0, nrow = length(dag.m[, 1]), ncol = length(dag.m[, 1]))
        diag(data.param.var) <- 1
        rownames(data.param.var) <- colnames(data.param.var) <- colnames(dag.m)
    }

    ## Due to historical reason this file is written using the transpose of the data.param matrix. Indeed, actual definition: children are on the rows and parents are on the column dag.m <- t(dag.m)
    ## data.param <- t(data.param)

    ##-----------------------------------------------------------------------------------------
    ## Creation of BUG file
    ##-----------------------------------------------------------------------------------------
    if (missing(bug.file)){
      bug.file <- "model.bug"
      newBugFile <- TRUE
    } else if (!file.exists(bug.file)) {
      newBugFile <- TRUE
    }
    if (newBugFile) {
        if (verbose == TRUE) {
            cat("Creation of the BUG file: '",bug.file,"'\n", sep = "")
        }

        ## erase existing file (if any)
        if (file.exists(bug.file)) {
            file.remove(bug.file)
        }

        ## create a new file
        sink(bug.file)
        cat("model{\n\n")

        ## start of the loop
        for (i in 1:length(dag.m[1, ])) {

            ## binomial nodes
            if (data.dists[i] == "binomial") {

                ## logistic regression
                cat("###-----------------------\n###Binomial nodes\n###-----------------------\n")
                cat(paste(names(dag.m[i, ])[i], " ~ dbern(p.", names(dag.m[i, ])[i], "); #Binary response\n", sep = ""))
                cat(paste("logit(", "p.", names(dag.m[i, ])[i], ") <- ", data.param[i, i], sep = ""))

                ## covariates
                for (j in 1:length(dag.m[1, ])) {
                  if (dag.m[i, j] == 1) {
                    cat(paste(" + ", data.param[i, j], "*", names(dag.m[i, ])[j], sep = ""))
                  }
                }

                cat(paste("; #Logistic regression \n\n"))

            }
            ## end of binomial nodes

            ## gaussian nodes
            if (data.dists[i] == "gaussian") {



                ## linear regression
                cat("###-----------------------\n###Gaussian nodes\n###-----------------------\n")
                cat(paste(names(dag.m[i, ])[i], " ~ dnorm(mu.", names(dag.m[i, ])[i], ", precision.", names(dag.m[i, ])[i], "); #Gausianne response\n", sep = ""))
                cat(paste("mu.", names(dag.m[i, ])[i], " <- ", data.param[i, i], sep = ""))

                ## covariates
                for (j in 1:length(dag.m[1, ])) {
                  if (dag.m[i, j] == 1) {
                    cat(paste(" + ", data.param[i, j], "*", names(dag.m[i, ])[j], sep = ""))
                  }
                }

                cat(paste("; #linear regression\n\n"))

                ## initialisation of the precision parameters
                cat(paste("precision.", names(dag.m[i, ])[i], " <- ", data.param.var[i, i], "; \n\n", sep = ""))


            }
            # EOF gaussian nodes

            ## poisson nodes
            if (data.dists[i] == "poisson") {

                ## poisson regression
                cat("###-----------------------\n###Poisson nodes\n###-----------------------\n")
                cat(paste(names(dag.m[i, ])[i], " ~ dpois(lambda.", names(dag.m[i, ])[i], "); #count response\n", sep = ""))
                cat(paste("log(", "lambda.", names(dag.m[i, ])[i], ") <- ", data.param[i, i], sep = ""))

                ## covariates
                for (j in 1:length(dag.m[1, ])) {
                  if (dag.m[i, j] == 1) {
                    cat(paste(" + ", data.param[i, j], "*", names(dag.m[i, ])[j], sep = ""))
                  }
                }

                cat(paste("; \n\n"))

            }

            ## multinomial nodes
            if (data.dists[i] == "multinomial") {

                l <- length(data.param.mult[i, ])

                ## logistic regression
                cat("###-----------------------\n###Multinomial nodes\n###-----------------------\n")
                cat(paste(names(dag.m[i, ])[i], "[1:", l, "] ~ dmulti(", names(dag.m[i, ])[i], ".p[],", sep = ""))

                cat(paste(names(dag.m[i, ])[i], ".N); #multinomial response\n", sep = ""))
                cat(paste(names(dag.m[i, ])[i], ".N <- ", data.param[i, i], ";\n", sep = ""))

                for (k in 1:l) {
                  cat(paste(names(dag.m[i, ])[i], ".p[", k, "]", " <- ", data.param.mult[i, k], ";\n", sep = ""))
                }

                cat(paste("\n\n"))

            }
        }

        ## close the model
        cat("}")
        ## closing file
        sink()
        ## closing file sink()

        if (verbose == TRUE) {
            cat("BUG '",bug.file,"' file created\n", sep = "")
        }
    }


    ##-----------------------------------------------------------------------------------------
    ## Simulation step
    ##-----------------------------------------------------------------------------------------

    if (simulate) {

        jagsverbose <- options(jags.pb = ifelse( verbose, "text", "none"))
        on.exit(jagsverbose)

        # Initial values
        init <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)

        jags <- jags.model(bug.file, inits = init, n.chains = n.chains, n.adapt = n.adapt, quiet = !verbose)

        res <- coda.samples(model = jags, variable.names = as.character(names(data.dists)), thin = n.thin, n.iter = n.iter)

        res <- do.call(rbind.data.frame, res)

        ## formatting
        for (i in 1:length(data.dists)) {
            if (data.dists[i] == "binomial") {
                res[, i] <- as.factor(res[, i])
            }
        }

        if(!keep.file) {
            if (file.exists(bug.file)) {
                file.remove(bug.file)
                }}

        return(as.data.frame(res))

    }
}  #EOF
