## plot-abn.R --- Author : Gilles Kratzer Last Modified on: 06/12/2016 Last Modified on: 10/03/2017 Last modification: 19.05.2017 Node color list Last mod: 13.06.2017 Arc direction Last mod:
## 18/07/2017

##-------------------------------------------------------------------------
## Function that plot a dag from a matrix or a formula
##-------------------------------------------------------------------------

plotabn <- function(dag.m = NULL, data.dists = NULL, markov.blanket.node = NULL, fitted.values.abn = NULL, 
                    fitted.values.abn.mle = NULL, digit.precision = 2, arc.strength = NULL, edgemode = "directed", 
    edgedir = "pc", node.fillcolor = "lightblue", edge.color = "black", edge.arrowwise = 0.5, fontsize.node = 10, 
    fontsize.edge = 5, plot = TRUE, node.fillcolor.list = NULL) {
    
    #while (!is.null(dev.list()))  dev.off()
    
    if(!is.null(markov.blanket.node) & ("multinomial" %in% (data.dists))) warning("Multinomial nodes are excluded from markov blanket computation.")
    
    ## for compatibility purpose
    dag <- dag.m
    if(inherits(x = dag,what = "abnLearned")){data.dists <- dag$score.cache$data.dists; dag <- dag$dag}
    
    group.var <- NULL
    name <- names(data.dists)
    
    
    ## dag transformation
    if (!is.null(dag)) {
        if (is.matrix(dag)) {
            ## run a series of checks on the DAG passed
            dag <- abs(dag)
            ## consistency checks
            diag(dag) <- 0
            dag[dag > 0] <- 1
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
    
    if (!plot) {
        dag[dag != 0] <- 1
        return(dag)
    }
    
    if (plot) {
        if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
            stop("library Rgraphviz is not available!\nRgraphviz is available as part of the bioconductor project - see http://www.bioconductor.org/install.\nRgraphviz is required when 'plot=TRUE'.")
        }
        
        if (edgemode == "undirected") {
            dag = dag + t(dag)
            dag[dag != 0] <- 1
        }
        
        
        if (edgedir == "pc") {
            dag <- t(dag)
        }
        if (edgedir == "cp") {
            dag <- (dag)
        }
        
        ## create an object graph
        am.graph <- new(Class = "graphAM", adjMat = (dag), edgemode = edgemode)
        
        ## close device dev.off()
        
        ## ======================================================================= Plotting options 1 Nodes options 2 Edges options =======================================================================
        
        ## ========= SHAPE =========
        
        ## Shape: plot differentially depending on the distribution
        
        shape <- list()
        
        for (i in 1:length(colnames(dag))) {
            if (data.dists[[i]] == "binomial") {
                shape <- cbind(shape, "box")
            }
            if (data.dists[[i]] == "gaussian") {
                shape <- cbind(shape, "circle")
            }
            if (data.dists[[i]] == "poisson") {
                shape <- cbind(shape, "ellipse")
            }
            if (data.dists[[i]] == "multinomial") {
                shape <- cbind(shape, "plaintext")
            }
        }
        
        class(shape) <- "character"
        names(shape) <- names(data.dists)
        
        ## ================= NODE FILLED COLOR =================
        
        if (exists("node.fillcolor")) {
            node.fillcolor <- rep(node.fillcolor, length(data.dists))
            class(node.fillcolor) <- "character"
            names(node.fillcolor) <- names(data.dists)
        }
        
        ## =============== MARKOV BLANKET ===============
        
        ## Markov Blanket: plot the MB of a given node
        
        if (!is.null(markov.blanket.node)) {
            
            ## Markov Blanket
            markov.blanket <- mb(dag, node = markov.blanket.node, data.dists = data.dists)
            
            ## Lists
            node.fillcolor <- list()
            names.node.fillcolor <- list()
            
            for (i in 1:length(colnames(dag))) {
                
                if (names(data.dists)[i] %in% as.list(markov.blanket.node)) {
                  
                  node.fillcolor <- cbind(node.fillcolor, "lightblue")
                  names.node.fillcolor <- cbind(names.node.fillcolor, names(data.dists)[i])
                  
                } else if (names(data.dists)[i] %in% as.list(unlist(markov.blanket))) {
                  
                  node.fillcolor <- cbind(node.fillcolor, "red")
                  names.node.fillcolor <- cbind(names.node.fillcolor, names(data.dists)[i])
                  
                } else {
                  
                  node.fillcolor <- cbind(node.fillcolor, "lightblue")
                  names.node.fillcolor <- cbind(names.node.fillcolor, names(data.dists)[i])
                }
            }
            
            
            class(node.fillcolor) <- "character"
            names(node.fillcolor) <- names.node.fillcolor
            
        }
        
        ## =============== Node color ===============
        
        if (!is.null(node.fillcolor.list)) {
            
            ## Lists
            node.fillcolor <- list()
            names.node.fillcolor <- list()
            
            for (i in 1:length(colnames(dag))) {
                
                if (names(data.dists)[i] %in% as.list(node.fillcolor.list)) {
                  
                  node.fillcolor <- cbind(node.fillcolor, "brown3")
                  names.node.fillcolor <- cbind(names.node.fillcolor, names(data.dists)[i])
                  
                } else {
                  
                  node.fillcolor <- cbind(node.fillcolor, "chartreuse3")
                  names.node.fillcolor <- cbind(names.node.fillcolor, names(data.dists)[i])
                }
            }
            
            class(node.fillcolor) <- "character"
            names(node.fillcolor) <- names.node.fillcolor
            
        }
        
        
        ## Edges options
        
        ## Edges names
        names.edges <- names(Rgraphviz::buildEdgeList(am.graph))
        
        ## =============== Fitted values ===============
        
        ## Plot the fitted values in abn as edges label
        
        if (!is.null(fitted.values.abn)) {
            tmp <- list()
            space <- "      "
            for (i in 1:length(fitted.values.abn)) {
                if (length(fitted.values.abn[[i]]) > 1) {
                  
                  if (data.dists[names(fitted.values.abn)[i]] == "gaussian" & length(fitted.values.abn[[i]]) > 2) {
                    
                    
                    for (j in 1:(length(fitted.values.abn[[i]]) - 2)) {
                      tmp <- cbind(tmp, paste(space, signif(fitted.values.abn[[i]][j + 1], digits = digit.precision)))
                    }
                    
                  } else {
                    
                    if (data.dists[names(fitted.values.abn)[i]] != "gaussian") {
                      for (j in 1:(length(fitted.values.abn[[i]]) - 1)) {
                        tmp <- cbind(tmp, paste(space, signif(fitted.values.abn[[i]][j + 1], digits = digit.precision)))
                      }
                    }
                  }
                }
            }
            
            label.edge <- tmp
            # return(label.edge)
            
        } else {
            
            label.edge <- rep(" ", length(names.edges))
            
        }
        
        ## ================= Fitted values MLE =================
        
        ## Plot the fitted values in abn as edges label (MLE output)
        
        if (!is.null(fitted.values.abn.mle)) {
            tmp <- list()
            space <- "      "
            for (i in 1:length(fitted.values.abn.mle)) {
                if (length(fitted.values.abn.mle[[i]]) > 1) {
                  
                  for (j in 1:(length(fitted.values.abn.mle[[i]]) - 1)) {
                    tmp <- cbind(tmp, paste(space, signif(fitted.values.abn.mle[[i]][j + 1], digits = digit.precision)))
                  }
                }
            }
            
            label.edge <- tmp
            # return(label.edge)
            
        } else {
            
            label.edge <- rep(" ", length(names.edges))
            
        }
        
        
        class(label.edge) <- "character"
        names(label.edge) <- names.edges
        
        ## =================== Arc Strength ===================
        
        ## Arc strength: plot the AS of the dag arcs
       
        if (is.matrix(arc.strength)) {
            if (edgemode != "undirected") {
            arc.strength <- t(arc.strength)
            
            min.as <- min(arc.strength[arc.strength > 0])
            max.as <- max(arc.strength[arc.strength > 0])
            
            
            arc.strength.norm <- (1/(max.as - min.as)) * (arc.strength - min.as)
            arc.strength.norm[arc.strength.norm < 0] <- 0
            lwd.edge <- list()
            for (i in 1:length(dag[1, ])) {
                for (j in 1:length(dag[1, ])) {
                  if (dag[i, j] == 1) {
                    lwd.edge <- cbind(lwd.edge, round(10 * arc.strength.norm[i, j]) + 1)
                  }
                }
            }
            
            # return(lwd.edge)
        }} else {
            
            lwd.edge <- rep(1, length(names.edges))
        }
        class(lwd.edge) <- "character"
        names(lwd.edge) <- names.edges
        
        ## ====== Plot ====== dev.off()
        Rgraphviz::plot(am.graph, attrs = list(node = list(fontsize = fontsize.node, fixedsize = FALSE), edge = list(arrowsize = edge.arrowwise, color = edge.color, fontsize = fontsize.edge)), nodeAttrs = list(fillcolor = node.fillcolor, 
            shape = shape), edgeAttrs = list(label = label.edge, lwd = lwd.edge))
    }
    
}  #EOF
