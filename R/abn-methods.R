###############################################################################
## abn-methods.R ---
## Author : Gilles Kratzer & Reinhard Furrer
## Document created : 21/05/2019
###############################################################################


##-------------------------------------------------------------------------
## abnDag
##-------------------------------------------------------------------------

# print

print.abnDag <- function(x, digits = 3L, ...){
  print(x$dag, digits = digits)
  cat("Class 'abnDag'.\n")
  invisible(x)
}

# summary

summary.abnDag <- function(object, ...) {
  su <- infoDag(object$dag)
  return(su)
}


# plot

plot.abnDag <- function(x, new=TRUE, ...){
  if (new) dev.new()
  on.exit(dev.flush())

  # Rgraphviz:
  mygraph <- new("graphAM", adjMat = t(x$dag), edgemode = "directed")
  g <- Rgraphviz::plot(x = mygraph)
  invisible(g)
}


##-------------------------------------------------------------------------
## abnCache
##-------------------------------------------------------------------------

# print

print.abnCache <- function(x, digits = 3, ...){

  cat("Number of nodes in the network: ",max(x$children), ".\n\n", sep='')
  if(x$method=="bayes"){
    cat("Distribution of the marginal likelihood: \n")
    print(summary(x[["mlik"]]), digits=digits)
  }

  if(x$method=="mle"){
    cat(" Distribution of the aic: \n")
    print(summary(x[["aic"]]), digits=digits)

    cat("\n Distribution of the bic: \n")
    print(summary(x[["bic"]]), digits=digits)

    cat("\n Distribution of the mdl: \n")
    print(summary(x[["mdl"]]), digits=digits)
  }
  invisible(x)
}

##-------------------------------------------------------------------------
## abnHeuristic
##-------------------------------------------------------------------------

# print

print.abnHeuristic <- function(x, digits = 2L, ...){
  cat("Best DAG' score found with",x$algo,"algorithm with", x$num.searches,"different searches limited to" , x$max.steps,"steps:\n")
  print(max(unlist(x$scores)), digits=digits)

  cat("\n Score distribution: \n")
  print(summary(unlist(x[["scores"]])), digits=digits)

  invisible(x)
}

# plot

plot.abnHeuristic <- function(x, ...){

  df <- unlist(x$scores)

  par(mfrow=c(1,2))
  plot(NULL, lty=1, xlab="Index of heuristic search", ylab="BN score", ylim = range(df), xlim = c(1,length(df)))
  for(i in 1:length(df)){
    if(sum(i==order(df, decreasing = FALSE)[1:10])){
      points(x=i,y=df[i], type="p", pch=19, col=rgb(0,0,1, 0.8),lwd = 2)
    } else {
      points(x=i,y=df[i], type="p", pch=19, col=rgb(0,0,0, 0.3))
    }
  }
  points(x = which.max(df), y = df[which.max(df)], col="red", pch=19)
  title("Networks final score")


  L <- (x$detailed.score)

  test <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))

  plot(NULL,lty=1, xlab="Number of Steps",ylab="BN score", ylim = range(test), xlim = c(1,length(test[,,1])))
  for(i in 1:length(L)){
    if(sum(i==order(df,decreasing = FALSE)[1:10])){
      points(x=1:(length(test[,,1])),y=test[1,,i], type="l", lty=1, col=rgb(0,0,1, 0.8),lwd = 2)
    } else {
      points(x=1:(length(test[,,1])),y=test[1,,i], type="l", lty=1, col=rgb(0,0,0, 0.17))
    }
  }
  lines(x=1:(length(test[,,1])),y=test[1,,which.max(df)], type="l", col="red", lwd=3)
  title("Networks score trajectory")
  invisible(x)
}

##-------------------------------------------------------------------------
## abnHillClimber
##-------------------------------------------------------------------------

# print

print.abnHillClimber <- function(x, digits = 3L, ...){
  print(x$consensus, digits = digits)
  cat("Consensus DAG from 'search.hillclimber'  (class 'abnHillClimber').\n")
  invisible(x)
}

# plot

plot.abnHillClimber <- function(x, new=TRUE, ...){

  if (new) dev.new()
  on.exit(dev.flush())

  # Rgraphviz
  mygraph <- new("graphAM", adjMat = x$consensus, edgemode = "directed")
  g <- Rgraphviz::plot(x = mygraph)
  invisible(g)
}


##-------------------------------------------------------------------------
## abnMostprobable
##-------------------------------------------------------------------------

# print

print.abnMostprobable <- function(x, digits = 3L, ...){

  print(x$dag, digits = digits)
  cat("Consensus DAG from 'mostprobable', can be use with 'fitabn'.\n")
  invisible(x)
}

# summary

summary.abnMostprobable <- function(object, ...){
  cat("Optimal DAG from 'mostProbable':\n")
  print(object$dag)
  cat( paste0("Calculated on ", dim(object$score.cache$data.df)[1], " observations.\n"))
  cat( paste0("(Cache length ", length(object$score.cache$mlik), '.)\n'))
  invisible( object)
}



# plot

plot.abnMostprobable <- function(x, new=TRUE, ...){

  if (new) dev.new()
  on.exit(dev.flush())

  # Rgraphviz:
  mygraph <- new("graphAM", adjMat = t(x$dag), edgemode = "directed")
  g <- Rgraphviz::plot(x = mygraph)
  invisible(g)
}

##-------------------------------------------------------------------------
## abnFit
##-------------------------------------------------------------------------

# print

print.abnFit <- function(x, digits = 3L, ...){

  if(x$method=="mle"){
    cat("The ABN model was fitted using an mle approach. The estimated coefficients are:\n\n")
    print(x$coef, digits=digits)
    cat("Number of nodes in the network: ",length(x$coef), ".\n", sep='')
  }

  if(x$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n\n")
    print(x$modes, digits=digits)
    cat("Number of nodes in the network: ",length(x$modes), ".\n", sep='')
  }

  invisible(x)
}

# summary

summary.abnFit <- function(object, digits = 3L, ...){

  if(object$method=="mle"){
    cat("The ABN model was fitted using an mle approach. The estimated coefficients are:\n")
    print(object$coef, digits=3)

    cat("Number of nodes in the network: ",length(object$modes), ".\n", sep='')

    cat("The AIC network score per node is: \n")
    print(unlist(object[["aicnode"]]), digits=digits)

    cat("\n The BIC network score per node is: \n")
    print(unlist(object[["bicnode"]]), digits=digits)

    cat("\n The MDL network score per node is: \n")
    print(unlist(object[["mdlnode"]]), digits=digits)
  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n")
    print(object$modes, digits=digits)

   cat("Number of nodes in the network: ",length(object$modes), ".\n\n", sep='')

   cat("The network score per node is:\n")
   print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}

# coef

coef.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){
  if(object$method=="mle"){
    cat("The ABN model was fitted using an mle approach. The estimated coefficients are:\n")
    print(object$coef, digits=digits)
  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n")
    print(object$modes, digits=digits)
  }

  invisible(object)

}



AIC.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The AIC network score per node is: \n")
    print(unlist(object[["aicnode"]]), digits=digits)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. AIC does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)


}

BIC.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The BIC network score per node is: \n")
    print(unlist(object[["bicnode"]]), digits=3)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. BIC does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}


logLik.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The loglikelihood network score per node is: \n")
    print(unlist(object[["mliknode"]]), digits=3)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. Loglikelihood does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}



family.abnFit <- function(object, ...){

  cat("All link functions are canonical: \n
      gaussian node = identy, binomial node = logit, Poisson node = log and multinomial node = logit.\n\n")

  print(unlist(object$abnDag$data.dists))

  invisible(object)
}


nobs.abnFit <- function(object, ...){
  nrow(object$abnDag$data.df)
}

plot.abnFit <- function(x, which ="abnFit", ...){

  if (which != "abnFit") stop('Function type not implemented yet. Use which="abnFit"')

  if (hasArg(fitted.values)) {
        g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, ...)
  } else {
    if(x$method=="mle"){
    g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, fitted.values = x$coef, ...)
  } else {
    g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, fitted.values = x$modes, ...)
  }
  }
  invisible(g)
}



## EOF
