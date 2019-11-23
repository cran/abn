################################################################################ 
## searchHeuristic.R 
## --- Author : Gilles Kratzer 
## Document created: 02/07/2018 
## Update: 21/05/2019

##-------------------------------------------------------------------------
## Heuristic search procedure
##-------------------------------------------------------------------------


searchHeuristic <- function(score.cache = NULL, score = "mlik", num.searches = 1, max.steps = 100, seed = 42, verbose = FALSE, start.dag = NULL, algo = "hc", tabu.memory = 10, temperature = 0.9, 
    ...) {
    
    ## method abnCache
    
    if (!inherits(score.cache,"abnCache")) {
        stop("score.cache should be an object of class 'abnCache' ")
    }
    
    data.dists <- score.cache$data.dists
    if(is.vector(score.cache$max.parents)){
        max.parents <- max(score.cache$max.parents)
    } else {
    max.parents <- score.cache$max.parents
    }
    dag.retained <- score.cache$dag.retained
    dag.banned <- score.cache$dag.banned
    
    ## function
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    
    n.var <- length(data.dists)
    
    ## output
    out.dags <- NULL
    out.scores <- list()
    
    
    out.detailed <- NULL
    
    
    ## seeding
    set.seed(seed = seed)
    
    for (searchIndex in 1:num.searches) {
        if (verbose) 
            cat("processing search...", searchIndex, "\n")
        
        temperature.update <- 1
        
        ## Initializing matrix
        dag.tmp <- matrix(data = 0, nrow = n.var, ncol = n.var)
        
        ## start zero matrix
        if (!is.null(start.dag)) {
            colnames(dag.tmp) <- rownames(dag.tmp) <- sample(1:n.var)
        }
        
        ## start random matrix
        if (is.null(start.dag)) {
            start.dag <- "random"
        }
        if (start.dag == "random") {
            vec.tmp <- c(rep(1, max.parents), rep(0, 2 * n.var - max.parents))
            for (lines in 1:(n.var - 1)) {
                dag.tmp[1 + lines, 1:lines] <- sample(vec.tmp)[1:lines]
            }
            colnames(dag.tmp) <- rownames(dag.tmp) <- sample(1:n.var)
        }
        if (is.matrix(start.dag)) 
            dag.tmp <- start.dag
        
        ## score init dag
        
        score.init <- vector(mode = "numeric", length = n.var)
        
        if (score %in% c("bic", "aic", "mdl")) {
            sc <- cbind(score.cache$node.defn[, as.numeric(colnames(dag.tmp))], -score.cache[[score]])
        } else {
            sc <- cbind(score.cache$node.defn[, as.numeric(colnames(dag.tmp))], score.cache[[score]])
        }
        
        for (lines in 1:n.var) {
            
            sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[lines]), ]
            
            score.init[lines] <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[lines, ])))), n.var + 1])
            
        }
        
        # dag.init <- dag.tmp start hc
        steps <- 1
        score.tmp <- score.init
        dag <- dag.tmp
        
        ## hill climbing algorithm
        if (algo == "hc") {
            
            
                out <- NULL
            
            
            while (steps < max.steps) {
                
                dag.tmp <- dag
                
                # change dag
                y <- resample(2:n.var, 1)
                x <- resample(1:(y - 1), 1)
                
                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1
                  
                } else {
                  
                  dag.tmp[y, x] <- 0
                }
                
                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])
                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }
                steps <- steps + 1
                
                
                  out <- cbind(out, sum(score.tmp))
                
                
            }  #eow
            
        }
        
        ## tabu algorithm
        if (algo == "tabu") {
            
            
                out <- NULL
            
            
            memory <- matrix(data = 0, nrow = tabu.memory, ncol = 2)
            prob <- c(0, 0)
            
            while (steps < max.steps) {
                
                dag.tmp <- dag
                
                ## tabu step
                while (sum(apply(tail(memory, tabu.memory), 1, function(x) identical((x), as.numeric(prob)))) > 0) {
                  # change dag
                  y <- resample(2:n.var, 1)
                  x <- resample(1:(y - 1), 1)
                  
                  prob <- c(x, y)
                  
                }
                
                memory <- rbind(memory, prob)
                
                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1
                  
                } else {
                  
                  dag.tmp[y, x] <- 0
                }
                
                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])
                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }
                steps <- steps + 1
                
                
                  out <- cbind(out, sum(score.tmp))
                
                
            }  #eow
        }
        
        ## simulated annealing (threashold acceptance)
        if (algo == "sa") {
            
            
                out <- NULL
            
            
            temperature.update <- 1
            
            
            while (steps < max.steps) {
                
                dag.tmp <- dag
                
                # change dag
                y <- resample(2:n.var, 1)
                x <- resample(1:(y - 1), 1)
                
                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1
                  
                } else {
                  
                  dag.tmp[y, x] <- 0
                }
                
                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])
                
                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                  
                }
                ## threshold acceptance (Dueck and Scheuer 1990) Dueck, G. and Scheuer, T. 'Threshold Accepting: A General Purpose Optimization Algorithm Appearing Superior to Simulated Annealing.' J. Comp.
                ## Phys. 90, 161-175, 1990.
                if (((score.tmp[y] - score.test) > 0) && (abs(score.tmp[y] - score.test) < (temperature.update * abs(score.init)/n.var))) {
                  
                  if (verbose) {
                    print("Accepting negative move")
                  }
                  
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }
                
                steps <- steps + 1
                temperature.update <- temperature.update * temperature
                
                
                  out <- cbind(out, sum(score.tmp))
                
                
            }  #eow
            
        }
        
        ## output of one search
        dag <- dag[as.character(1:n.var), as.character(1:n.var)]
        colnames(dag) <- rownames(dag) <- names(data.dists)
        
        out.dags[[searchIndex]] <- (dag)
        out.scores[[searchIndex]] <- sum(score.tmp)
        
        
        out.detailed[[searchIndex]] <- out
        
        
    }  #eos
    
    
    ## return
    
    if (score %in% c("bic", "aic", "mdl")) {
            
            out <- list(dags = out.dags, 
                 scores = lapply(X = out.scores, FUN = function(x) {
                     -x
                 }), detailed.score = lapply(X = out.detailed, FUN = function(x) {
                     -x
                 }), score = score, score.cache = score.cache, num.searches = num.searches, max.steps = max.steps,algorithm = algo
            )
            
            class(out) <- c("abnHeuristic", "mle")
            
            return(out)
            
            #return(list(dag = out.dags[[which.max(unlist(out.scores))]], score = -max(unlist(out.scores))))
        
    } else {
        out <- list(dags = out.dags, scores = out.scores, detailed.score = out.detailed,
                    score = score, score.cache = score.cache, num.searches = num.searches, max.steps = max.steps,algorithm = algo)
        class(out) <- c("abnHeuristic", "bayes")
            return(out)

            #return(list(dag = out.dags[[which.max(unlist(out.scores))]], score = max(unlist(out.scores))))
        }

    
}  #eof
