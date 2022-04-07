##############################################################################
# abn-internal.R

abn.Version <- function(what=c('abn','system')) {
  what <- match.arg(what)
  if (what %in% 'system') {

    list(R=R.version.string,
         abn=substr(abn.version$version.string, 13, 32),
         gsl=ifelse(R.version$os=="linux-gnu", system('gsl-config --version', intern = TRUE), "NA (?)"),
         JAGS=rjags::jags.version(),
         INLA=ifelse(requireNamespace("INLA", quietly = TRUE),
           INLA::inla.version("version"), "not available")
    )

  } else {
    release <- utils::packageDescription("abn",field="Version")
    date <- utils::packageDescription("abn",field="Date")
    list(status="",
        major=sub("-","",substr(release,1,4)),
        minor=substr(sub("-","",substr(release,5,7)),1,1),
         year=substr(date,1,4),
         month=substr(sub("20..-","",date),1,2),
         day=sub("20..-..-","",date),
         version.string= paste("abn version ",
                             utils::packageDescription("abn",field="Version")," (",
                             utils::packageDescription("abn",field="Date"),")",sep="")
  )
  }
}

abn.version <- abn.Version()
class(abn.version) <- "simple.list"




".onAttach" <- function (lib, pkg) {
  packageStartupMessage(abn.version$version.string," is loaded.")
}

##-------------------------------------------------------------------------
## Internal function that call multiple times strsplit() and remove space
##-------------------------------------------------------------------------

strsplits <- function(x, splits, ...) {
    for (split in splits) {
        x <- unlist(strsplit(x, split, ...))
    }
    x <- gsub(" ", "", x, fixed = TRUE)  #remove space
    return(x[!x == ""])  # Remove empty values
}

##-------------------------------------------------------------------------
## Internal function that produce a square matrix length(name) with {0,1} depending on f. f have to start with ~ terms are entries of name terms are separated by + term1 | term2 indicates
## col(term1) row(term2) puts a 1 term1 | term2:term3: ... : is used as a sep . = all terms in name
##-------------------------------------------------------------------------

formula.abn <- function(f, name) {

    name_orignial <- name

    f <- as.character(f)

    ## tests for consistence ---------------------------------------------------------------------- start as a formula
    if (!grepl("~", f[1], fixed = T)) {
        stop("DAG specifications should start with a ~")
    }

    ## transformation name + or | or : or . or name to name_name
    if (sum((c("+", "|", ":", ".") %in% unlist(strsplit(name, split = c(""))))) != 0) {
        for (i in 1:length(name)) {
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("+")) != 0) {
                f[[2]] <- gsub(name[i], gsub("+", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("+", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c("|")) != 0) {
                f[[2]] <- gsub(name[i], gsub("|", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub("|", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(":")) != 0) {
                f[[2]] <- gsub(name[i], gsub(":", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(":", "_", name[i], fixed = TRUE)
            }
            if (sum(unlist(strsplit(name[i], split = c(""))) %in% c(".")) != 0) {
                f[[2]] <- gsub(name[i], gsub(".", "_", name[i], fixed = TRUE), f[[2]], fixed = TRUE)
                name[i] <- gsub(".", "_", name[i], fixed = TRUE)
            }
        }
    }

    ## collapse name
    name.c <- paste(name, collapse = ":")
    ## Split by terms
    f.p <- strsplit(x = f[[2]], split = "+", fixed = TRUE)

    ## nothing more than name variable in the dag formula
    tmp.test <- strsplits(x = f[[2]], splits = c("+", "|", ":", "."), fixed = TRUE)
    if (sum(!(tmp.test %in% name)) != 0) {
        stop("DAG formulation contains some variables not in provided names")
    }
    ## End of tests for consistence ----------------------------------------------------------------

    ## creat the void matrix
    out <- matrix(data = 0, nrow = length(name), ncol = length(name))

    ## delete all spaces
    f.p <- gsub(" ", "", f.p[[1]], fixed = TRUE)

    ## replace '.' by all names
    f.p.completed <- gsub(".", name.c, f.p, fixed = TRUE)

    ## atomization of left term


    ## contruction of the output matrix
    for (i in 1:length(f.p)) {
        tmp <- f.p.completed[i]

        ## forget unique terms -> test for |
        if (grepl("|", tmp, fixed = TRUE)) {

            ## split wrt |
            tmp.p <- strsplit(x = tmp, split = "|", fixed = TRUE)

            ## test for multiple terms and contruction of the list first term
            if (grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- strsplit(x = tmp.p[[1]][1], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][1])) {
                tmp.p.p.1 <- tmp.p[[1]][1]
            }

            ## test for multiple terms and contruction of the list second term
            if (grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- strsplit(x = tmp.p[[1]][2], split = ":", fixed = TRUE)
            }
            if (!grepl(":", tmp.p[[1]][2])) {
                tmp.p.p.2 <- tmp.p[[1]][2]
            }

            ## loop over the
            for (j in 1:length(tmp.p.p.1[[1]])) {
                for (k in 1:length(tmp.p.p.2[[1]])) {
                  ## update of matrix
                  out[grep(tmp.p.p.1[[1]][j], name), grep(tmp.p.p.2[[1]][k], name)] <- 1

                }
            }
        }

    }

    ## avoid auto dependance
    diag(out) <- 0

    ## only 0 and 1
    out[out > 1] <- 1

    ## naming
    colnames(out) <- name_orignial
    rownames(out) <- name_orignial
    ## output
    return(out)
}

#####################################################################################################
################### a set of simple commonsense validity checks on the data.df and data.dists arguments
check.valid.data <- function(data.df = NULL, data.dists = NULL, group.var = NULL) {

    ## check data is in a data.frame
    if (!is.data.frame(data.df)) {
        stop("The data must be in a data.frame")
    }

    ## check data for missing values
    if (sum(complete.cases(data.df)) != dim(data.df)[1]) {
        stop("The data contains missing values! These must be removed.")
    }

    ## check that distributions are in a list
    if (!is.list(data.dists)) {
        stop("data.dist must be a list")
    }

    if (!is.null(group.var)) {
        ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
        data.df <- data.df[, -which(names(data.df) == group.var)]
    }

    if (length(names(data.dists)) != length(names(data.df))) {
        stop("data.dists must have named entries")
    }

    ## check names in list are in correct order and exact match to data.frame
    for (i in 1:dim(data.df)[2]) {
        if (names(data.dists)[i] != names(data.df)[i]) {
            stop("names in list must match names in data.frame")
        }
    }

    ## check names of likelihood function are valid
    allowed.dists <- c("gaussian", "binomial", "poisson")
    ## n.b. these must match inla() family=''
    for (i in 1:length(data.dists)) {
        if (length(which(allowed.dists %in% data.dists[[i]])) != 1) {
            ## each variable must have one of the allowed.dists
            message <- paste("Each variable must have a distribution of either", paste(allowed.dists, collapse = ", "))
            cat(message, "\n")
            stop("")
        }
    }

    binomial.vars.indexes <- NULL
    poisson.vars.indexes <- NULL
    gaussian.vars.indexes <- NULL

    ## check that data is consistent with distribution given for each variable
    for (i in 1:dim(data.df)[2]) {
        cur.var <- data.df[, i]
        if (data.dists[[i]] == "gaussian") {
            if (is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must not be a factor.\n")
                stop("")
            }
            if (length(unique(cur.var)) <= 2) {
                cat((names(data.df)[i]), "is invalid as it has two or less unique values!\n")
                stop("")
            }
            gaussian.vars.indexes <- c(gaussian.vars.indexes, i)
        }
        if (data.dists[[i]] == "binomial") {
            if (!is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must be a factor\n")
                stop("")
            }
            if (length(unique(cur.var)) < 2) {
                cat((names(data.df)[i]), "is invalid as it must be binary with both cases being observed.\n")
                stop("")
            }
            if (length(unique(cur.var)) > 2) {
                cat((names(data.df)[i]), "is invalid as it must be binary. Multi-category variables should be split into separate binary variables.\n")
                stop("")
            }
            binomial.vars.indexes <- c(binomial.vars.indexes, i)
        }

        if (data.dists[[i]] == "poisson") {
            if (is.factor(cur.var)) {
                cat((names(data.df)[i]), "is invalid - it must not be a factor\n")
                stop("")
            }
            if (length(unique(cur.var)) <= 2) {
                cat((names(data.df)[i]), "is invalid as it has two or less unique values!")
                stop("")
            }
            poisson.vars.indexes <- c(poisson.vars.indexes, i)
        }
    }
    ## return the indexes of any binary variables
    return(list(gaus = gaussian.vars.indexes, bin = binomial.vars.indexes, pois = poisson.vars.indexes))

}  #end of check.valid.data()

######################################### a set of simple commonsense validity checks on the directed acyclic graph definition matrix
check.valid.dag <- function(dag = NULL, data.df = NULL, is.ban.matrix = FALSE, group.var = NULL) {

    if (!is.null(group.var)) {
        ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY

        if (is.null(data.df)) stop("When specifying 'group.var', 'data.df' argument is required as well.")
        data.df <- data.df[, -which(names(data.df) == group.var)]
    }


    ## if dag null then create unlimited - empty - network want ban matrix
    if (is.null(dag)) {
        dag <- matrix(rep(0, dim(data.df)[2]^2), ncol = dim(data.df)[2])
        ## names must be set
        colnames(dag) <- rownames(dag) <- names(data.df)
        return(dag)
    }

    if (!is.matrix(dag)) {
        if (grepl("~", as.character(dag)[1], fixed = T)) {
            dag <- formula.abn(f = dag, name = names(data.df))
            return(dag)
        } else {
            stop("'dag' specification must either be a matrix or a formula expresion")
        }
    } else {
        if (dim(dag)[1] != dim(dag)[2])  stop("Matrix 'dag' is not square.")
        }


    ## check data for missing names
    if (is.null(colnames(dag)) || is.null(rownames(dag))) {
        if (!is.null(data.df)) {
            if (dim(dag)[1] != dim(data.df)[2]) {  stop("'dag' as dimension inconsistent with columns of 'data.df'")
#                print(dag)
#                print(data.df)
                }
           colnames(dag) <- rownames(dag) <- names(data.df)
        } else {
            stop("'dag' must have both row and column names set or a named dataset has to be provided")
        }
    }
    ## check dimension
    if (!is.null(data.df)) {
      if (dim(dag)[1] != dim(data.df)[2] || dim(dag)[2] != dim(data.df)[2]) {
         stop("'dag' as dimension inconsistent with data.df")
      }
    }

    ## check binary
    for (i in 1:dim(dag)[1]) {
       for (j in 1:dim(dag)[2]) {
          if ( abs(dag[i, j]) > 1e-8 && abs(dag[i, j]-1) > 1e-8)      stop("'dag' must comprise only 1's or 0's")
       }
    }
    # if (any( c( (dag != 0) && (dag != 1))))                stop("'dag' must comprise only 1's or 0's")


    ## check diagnonal and cycles - but ignore these checks for a ban matrices
    if (!is.ban.matrix) {

        if (any( diag(dag) != 0))                 stop("'dag' is not a valid DAG - a child cannot be its own parent!")



        ## coerce to int for sending to C$ number of cols (or rows)
        ## this creates one long vector - filled by cols from dag = same as internal C reprentation so fine.
        res <- .Call("checkforcycles", as.integer(dag),  dim(dag)[1], PACKAGE = "abn")
        if (res!=0) stop("'dag' contains at least one cycle.")
    }

    return( dag)


}


######################################### a set of simple checks on the list given as parent limits
check.valid.parents <- function(data.df = NULL, max.parents = NULL, group.var) {
    ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    if (!is.null(group.var)) {
        data.df <- data.df[, -which(names(data.df) == group.var)]
    }
    # print(data.df);print(max.parents); if a constant then make integer vector
    if (is.numeric(max.parents) && length(max.parents) == 1) {
        return(as.integer(rep(max.parents, dim(data.df)[2])))
    }

    ## if a list must be named list with names as in original data
    if (is.list(max.parents) && length(max.parents) == dim(data.df)[2]) {
        for (i in 1:dim(data.df)[2]) {
            if (names(max.parents)[i] != names(data.df)[i]) {
                stop("names in max.parents list must match names in data.frame data.df")
            }
        }
        if (!is.numeric(unlist(max.parents))) {
            stop("max.parents is not valid - must be numeric")
        }
        max.parents.int <- unlist(max.parents)
        if (length(max.parents.int) != dim(data.df)[2]) {
            stop("max.parents list is wrong length")
        }
        max.parents.int <- as.integer(max.parents.int)
        return(max.parents.int)
    }

    stop("'max.parents' is not valid: length data: ",dim(data.df)[2],
         ", length max.parents: ",length(max.parents))

}

######################################### a set of simple checks on the list given as parent limits
check.which.valid.nodes <- function(data.df = NULL, which.nodes = NULL, group.var) {
    ## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    if (!is.null(group.var)) {
        data.df <- data.df[, -which(names(data.df) == group.var)]
    }
    ## if null then assume ALL nodes
    if (is.null(which.nodes)) {
        which.nodes <- 1:dim(data.df)[2]
        return(as.integer(which.nodes))
    }

    if (is.numeric(which.nodes) && max(which.nodes) <= dim(data.df)[2] && min(which.nodes) >= 1 && length(which.nodes) <= dim(data.df)[2]) {
        return(as.integer(which.nodes))
    } else {
        stop("which.nodes is invalid")
    }

}

######################################### a simple check on the grouping variable
check.valid.groups <- function(group.var, data.df, cor.vars) {

    if (is.null(group.var))
        {
            return(list(data.df = data.df, grouped.vars = as.integer(c(-1)), group.ids = as.integer(rep(-1, dim(data.df)[1]))))
        }  ## have no groups so just return dummy values
    if (!(is.character(group.var) && (length(group.var) == 1))) {
        stop("name of group variable is not a character?!")
    }
    if (!length(which(group.var %in% names(data.df) == TRUE))) {
        stop("name of group variable does not match any of those in data.df")
    }
    ## get group id data
    group.var.vals <- data.df[, group.var]
    ## drop the group variable from original data.frame and overwrite
    data.df <- data.df[, -which(names(data.df) == group.var)]

    ## have groups so some checks

    if (is.factor(group.var.vals) && length(group.var.vals) == dim(data.df)[1] && length(unique(group.var.vals)) > 1) {
        ## is factor and of correct length and at least two groups
    } else {
        stop("grouping variable must be: i) a factor; ii) same length as data.df; and iii) contain more than one group")
    }

    ## get group memberships in terms of ints
    group.var <- as.integer(group.var.vals)

    ## now find out which variables the grouping is to be applied to
    var.noms <- names(data.df)
    if (length(which(cor.vars %in% var.noms == TRUE)) != length(cor.vars)) {
        stop("variables in cor.vars do not match those in data.df")
    }

    if (max(table(cor.vars)) > 1) {
        stop("have repeated variables in cor.vars!")
    }

    ## to get to here group.var must be ok and also variable names so return integer code for the variables get the index in names(data.df) for each variable and then sort into order
    cor.var.indexes <- as.integer(sort(match(cor.vars, var.noms)))


    return(list(data.df = data.df, grouped.vars = cor.var.indexes, group.ids = group.var))

}

########################################## create ordered vector with integers denoting the distribution
get.var.types <- function(data.dists = NULL) {
    store <- rep(NA, length(data.dists))

    for (i in 1:length(data.dists)) {
        if (data.dists[[i]] == "binomial") {
            store[i] <- 1
        }
        if (data.dists[[i]] == "gaussian") {
            store[i] <- 2
        }
        if (data.dists[[i]] == "poisson") {
            store[i] <- 3
        }
    }

    return(store)

}
########################################## tidy up cache
tidy.cache <- function(thecache) {
    if (!is.null(thecache[["error.indexes"]])) {
        error.combs <- thecache[["error.indexes"]]
        corrected <- list()
        corrected[["children"]] <- thecache[["children"]][-error.combs]
        corrected[["node.defn"]] <- thecache[["node.defn"]][-error.combs, ]
        corrected[["mlik"]] <- thecache[["mlik"]][-error.combs]
        ## return new cache with appropriate nodes removed
        return(corrected)
    }

}

