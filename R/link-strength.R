################################################################################ link-strength.R --- Author : Gilles Kratzer Document created: 09/05/2017 -19/06/2017: ls ls.pc

##-------------------------------------------------------------------------
## Function that return estimation of link strength for a Bayesian Network
##-------------------------------------------------------------------------

linkStrength <- function(dag, data.df = NULL, data.dists = NULL, method = c("mi.raw", "mi.raw.pc", "mi.corr", "ls", "ls.pc", "stat.dist"), discretization.method = "doane") {

    group.var <- NULL
    if (is.matrix(dag)) {
        name <- names(dag)
    } else {
        name <- names(data.dists)
    }
    # Contains Rgraphviz stuff

    ## dag transformation
    if (!is.null(dag)) {
        if (is.matrix(dag)) {
            ## run a series of checks on the DAG passed
            dag <- abs(dag)
            diag(dag) <- 0
            dag <- check.valid.dag(dag.m = dag, is.ban.matrix = FALSE, group.var = group.var)
            ## naming
            if (is.null(colnames(dag))) {
                colnames(dag) <- name
                rownames(dag) <- name
            }
        } else {
            if (grepl("~", as.character(dag)[1], fixed = T)) {
                dag <- formula.abn(f = dag, name = name)
                ## run a series of checks on the DAG passed
                dag <- check.valid.dag(dag.m = dag, is.ban.matrix = FALSE, group.var = group.var)
            }
        }
    } else {
        stop("Dag specification must either be a matrix or a formula expression")
    }


    ## rows and columns
    n.row <- length(dag[1, ])
    n.col <- length(dag[, 1])


    ## MUTUAL INFORMATION

    if (method == "mi.raw" | method == "mi.raw.pc")
        {

            mi.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    mi.m[j, i] <- miData(freqs.table = discretization(data.df = data.df[, c(i, j)], data.dists = data.dists[c(i, j)], discretization.method), method = method)

                  }
                }
            }

            return(mi.m)
        }  #EOF: mi.raw

    if (method == "mi.corr")
        {

            mi.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    mi.m[j, i] <- -0.5 * log(1 - cor(x = data.df[, j], y = data.df[, i])^2)
                  }
                }
            }

            return(mi.m)
        }  #EOF: mi.theo

    if (method == "ls")
        {
            ## Formula: I(X;Y|Z)=H(X,Z)+H(Y,Z)-H(X,Y,Z)-H(Z)=H(X|Z)-H(X|Y,Z) with Z all other parents of Y i -> j = p -> c
            ls.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {
                    parent.list <- t(dag[j, ])
                    parent.list[parent.list != 0] <- 1

                    parent.list[j] <- 0
                    parent.list[i] <- 0

                    parent.list.x <- parent.list
                    parent.list.x[i] <- 1
                    parent.list.y <- parent.list
                    parent.list.y[j] <- 1

                    if (sum(parent.list) != 0) {
                      # ls.m <- rbind(parent.list)
                      ls.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method))
                    } else {
                      ls.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method))
                    }
                  }
                }
            }
            return(ls.m)
        }  #EOF: link strength


    if (method == "ls.pc")
        {
            ## Formula: I(X;Y|Z)/H(Y|Z)=(H(X,Z)+H(Y,Z)-H(X,Y,Z)-H(Z))/H(Y|Z)=(H(X|Z)-H(X|Y,Z))/H(Y|Z). H(Y|Z)=H(Y,Z)-H(Z)

            ls.pc.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {
                    parent.list <- t(dag[j, ])
                    parent.list[parent.list != 0] <- 1

                    parent.list[j] <- 0
                    parent.list[i] <- 0

                    parent.list.x <- parent.list
                    parent.list.x[i] <- 1
                    parent.list.y <- parent.list
                    parent.list.y[j] <- 1
                    parent.list.xy <- parent.list.y
                    parent.list.xy[i] <- 1

                    if (sum(parent.list) != 0) {
                      # ls.m <- rbind(parent.list)
                      ls.pc.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.xy)], data.dists = data.dists[as.logical(parent.list.xy)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method))
                      tmp <- (entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method)))
                      if (tmp != 0) {
                        ls.pc.m[j, i] <- ls.pc.m[j, i]/tmp
                      }

                    } else {
                      ls.pc.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.xy)], data.dists = data.dists[as.logical(parent.list.xy)], discretization.method = discretization.method))
                      tmp <- (entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)))
                      if (tmp != 0) {
                        ls.pc.m[j, i] <- ls.pc.m[j, i]/tmp
                      }
                    }
                  }
                }
            }
            return(ls.pc.m)
        }  #EOF: link strength


    if (method == "stat.dist")
        {
            ## Formula: statistical distance = 1- MI(X,Y)/max(H(X), H(Y))

            stat.dist <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    stat.dist[j, i] <- miData(freqs.table = discretization(data.df = data.df[, c(i, j)], data.dists = data.dists[c(i, j)], discretization.method), method = "mi.raw")

                  }
                }
            }
            return(stat.dist)
        }  #EOF: statistical distance


    # stat.dist <- matrix(data = 0,nrow = n.row,ncol = n.col) for(i in 1:n.col){ for(j in 1:n.row){ if(dag[j,i]!=0){ stat.dist[j,i] <- mi.data(freqs.table = discretization(data.df =
    # data.df[,c(i,j)],data.dists = data.dists[c(i,j)],discretization.method), method = method) #stat.dist[j,i] <- stat.dist[j,i]/max(entropy.data(freqs.table = discretization(data.df =
    # data.df[,i],data.dists=data.dists[i],discretization.method = discretization.method)), entropy.data(freqs.table = discretization(data.df =
    # data.df[,j],data.dists=data.dists[j],discretization.method = discretization.method))) #stat.dist[j,i] <- 1-stat.dist[j,i] }}} return(mi.m) }#EOF: statistical distance


}  #EOF
