## plot-abn.R --- Author : Gilles Kratzer Last Modified on: 06/12/2016 Last Modified on: 10/03/2017 Last modification: 19.05.2017 Node color list Last mod: 13.06.2017 Arc direction Last mod:
## 18/07/2017
## major rewrite rf 2021-04

# for final submission elimiate and # `print()` lines


plotabn <- function(...) {
    .Deprecated("plotAbn", msg="'plotabn' is deprecated.\n Use 'plotAbn' instead but note that arguments have slightly changed.")
    dots <- list(...)
    if (!is.null(dots$dag.m)) {
        dag <- dots$dag.m
        dots$dag.m <- NULL
        do.call('plotAbn', c(dag, dots))
    } else  plotAbn(...)
}


plotAbn <- function(dag, data.dists=NULL, markov.blanket.node=NULL,
                    fitted.values=NULL, digits=2, edge.strength=NULL, edge.strength.lwd=5,
                    edge.direction="pc", edge.color="black", edge.linetype="solid", edge.arrowsize=0.6, edge.fontsize=node.fontsize,
                    node.fontsize=12, node.fillcolor=c("lightblue","brown3","chartreuse3"),
                    node.fillcolor.list=NULL, node.shape=c("circle","box","ellipse","diamond"),
                    plot=TRUE , ... )       {

    # Actually, the plot argument is wrong! i do not need the adjacency structure only. I need all but the plotting. i.e., all but the rendering of the graph.


    # The following is not relevant. The nodes are calculated via mb. They are not colored.
    #    if(!is.null(markov.blanket.node) & ("multinomial" %in% (data.dists))) warning("Multinomial nodes are excluded from Markov blanket computation.")

    ## for compatibility purpose
    if(inherits(x=dag, what="abnLearned")){
        data.dists <- dag$score.cache$data.dists;
        dag <- dag$dag
    }
  if(inherits(x=dag, what="abnFit")){
    data.dists <- dag$abnDag$data.dists
    dag <- dag$abnDag$dag
  }
  if (is.null(data.dists)) stop("'data.dist' need to be provided.")
    name <- names(data.dists)


    ## dag transformation
    if (!is.null(dag)) {
        if (is.matrix(dag)) {
            ## run a series of checks on the DAG passed
            dag <- abs(dag)
            ## consistency checks
            diag(dag) <- 0
            dag[dag > 0] <- 1

            ## naming
            if (is.null(rownames(dag))) {
                colnames(dag) <- name
                rownames(dag) <- name
            }
            dag <- check.valid.dag(dag=dag, is.ban.matrix=FALSE, group.var=NULL)
        } else {
            if (grepl("~", as.character(dag)[1], fixed=T)) {
                dag <- formula.abn(f=dag, name=name)
                ## run a series of checks on the DAG passed
                dag <- check.valid.dag(dag=dag, is.ban.matrix=FALSE, group.var=NULL)
            }
        }
    } else {
        stop("'dag' specification must either be a matrix or a formula expression")
    }

    # contains Rgraphviz
    if (edge.direction == "undirected") {
        dag=dag + t(dag)
        dag[dag != 0] <- 1     # this should not be necessary!
    }


    ## create an object graph
    am.graph <- new(Class="graphAM", adjMat=dag,
                    edgemode=ifelse(edge.direction=="undirected","undirected","directed"))

    ## ========= SHAPE =========
    ## Shape: plot differential depending on the distribution
    node.shape <- rep(node.shape, 4)
    shape <- rep(node.shape[1], length(data.dists) )
    shape[data.dists == "binomial"] <- node.shape[2]
    shape[data.dists == "poisson"] <- node.shape[3]
    shape[data.dists == "multinomial"] <- node.shape[4]
    names(shape) <- names(data.dists)

    ## ================= NODE FILLED COLOR =================
    ## fill with default value, change if MB or fillcolor.list is requested
    fillcolor <- rep(node.fillcolor[1], length(data.dists))
    names(fillcolor) <- names(data.dists)

    ## =============== MARKOV BLANKET ===============
    ## Markov Blanket: plot the MB of a given node
    if (!is.null(markov.blanket.node)) {
        markov.blanket <- mb( dag, node=markov.blanket.node, data.dists=data.dists)
        fillcolor[ names(data.dists) %in%  markov.blanket]  <- node.fillcolor[3]
        fillcolor[ names(data.dists) %in%  markov.blanket.node]  <- node.fillcolor[2]

    } else    if (!is.null(node.fillcolor.list)) {
        fillcolor[ names(data.dists) %in%  node.fillcolor.list] <- node.fillcolor[2]
    }

    names.edges <- names(Rgraphviz::buildEdgeList(am.graph))

    ## =============== Fitted values ===============
    ## Plot the fitted values in abn as edges label
#    print(names.edges)
    if (!is.null(fitted.values)) {
        space <- "   "
        edge.label <- c()
        for (i in 1:length(fitted.values)) {
            if ((length(fitted.values[[i]]) > 1)& (data.dists[names(fitted.values)[i]] != "gaussian")) {
                for (j in 1:(length(fitted.values[[i]]) - 1))
                    edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][j + 1], digits=digits)))
            } else if ((length(fitted.values[[i]]) > 2)& (data.dists[names(fitted.values)[i]] == "gaussian")){
                for (j in 1:(length(fitted.values[[i]]) - 2))
                    edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][j + 1], digits=digits)))
            }
        }
    } else  edge.label <- rep(" ", length(names.edges))
    names(edge.label) <- names.edges


    ## =================== Arc Strength ===================
    ## Arc strength: plot the AS of the dag arcs
#    if (is.matrix(edge.strength) & (edge.direction != "undirected")) {
    if (is.matrix(edge.strength)) {
        if (any(edge.strength<0)) stop("'edge.strength' should be positive")
        if (any(edge.strength[dag ==0] >0)) stop("'edge.strength' does not match dag")
        min.as <- min(edge.strength[edge.strength > 0])
        max.as <- max(edge.strength[edge.strength > 0])

        edge.strength.norm <- (edge.strength - min.as)/(max.as - min.as)
        edge.strength.norm[edge.strength.norm < 0] <- 0
        edge.lwd <- list()
        for (i in 1:length(dag[1, ])) {
            for (j in 1:length(dag[1, ])) {
                if (dag[i, j] == 1) {
                    edge.lwd <- cbind(edge.lwd, round(edge.strength.lwd * edge.strength.norm[i, j]) + 1)
                }
            }
        }
    } else {
        edge.lwd <- rep(1, length(names.edges))
    }
    class(edge.lwd) <- "character"
    names(edge.lwd) <- names.edges

    ## ====== Plot ======
    attrs <- list(graph=list(rankdir="BT"),
                  node=list(fontsize=node.fontsize, fixedsize=FALSE),
                  edge=list(arrowsize=edge.arrowsize, color=edge.color, lty=edge.linetype, fontsize=edge.fontsize))
    nodeAttrs <- list(fillcolor=fillcolor, shape=shape)
    edgeAttrs <- list(label=edge.label, lwd=edge.lwd)
#     print(edgeAttrs)
#    if (all(shape %in% c("circle","box","ellipse")))  {
    am.graph <- layoutGraph(am.graph, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)

    if (edge.direction == "pc")  {     # specify appropriate direction!
      edgeRenderInfo(am.graph) <- list(arrowtail="open")
      edgeRenderInfo(am.graph) <- list(arrowhead="none")
#      edgeRenderInfo(am.graph) <- list(direction=NULL)# MESSES up!!! not needed.
    }
    edgeRenderInfo(am.graph) <- list(lwd=edge.lwd)

  #  if (plot) renderGraph(am.graph, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
    if (plot) renderGraph(am.graph, ...)

#   } else {
#        am.graph <- layoutGraph(am.graph, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs, ...)
        # the following does not work in R
#        edgeRenderInfo(am.graph)[["direction"]] <- "back"
        # hence
#        warning("edge.direction='pc' is not working with diamond shapes.")
#        edgeRenderInfo(am.graph) <- list(lwd=edge.lwd)
#        if (plot) renderGraph(am.graph,attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
#    }

    invisible(am.graph)
}  #EOF
