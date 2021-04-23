###############################################################################
## fitabn.R --- 
## Author          : Gilles Kratzer
## Document created: 13/02/2017
## Last modified   : 28/02/2017 High lvl optimizations
## Last modified   : 06/03/2017 Pvalue
## Last modified   : 20/10/2017 Multinomial 
## Last modified   : 23/10/2017 Multinomial naming if clauses
###############################################################################

## fit a given DAG to data
fitAbn.mle <- function(dag = NULL,
                       data.df = NULL, 
                       data.dists = NULL,
                       adj.vars = NULL,
                       cor.vars = NULL,
                       centre = TRUE,
                       maxit = 100, 
                       tol = 10^-11,
                       verbose = FALSE, seed = 9062019){

  set.seed(seed)
  
  n <- length(data.dists)
  nobs <- dim(data.df)[1]
  
  #test for dag
  if(!is.null(dag)){
    if(is.matrix(dag)){
      ## run a series of checks on the DAG passed
      dag <- check.valid.dag(dag.m=dag,data.df=data.df,is.ban.matrix=FALSE,group.var=NULL)
    } else {
      if(grepl('~',as.character(dag)[1],fixed = TRUE)){
        dag <- formula.abn(f = dag,name = colnames(data.df))
        ## run a series of checks on the DAG passed
        dag <- check.valid.dag(dag.m=dag,data.df=data.df,is.ban.matrix=FALSE,group.var=NULL)
      }
    }}
  else {
    stop("Dag specification must either be a matrix or a formula expression")
  }
  
  
  ##some bunch of test
  
  #test same order for data.frame and data.dist
  if(is.null(names(data.dists)))stop("data.dist is not a named vector")
  if(is.null(names(data.df)))stop("data.df is not a named data frame")
  if(is.null(colnames(dag)) | is.null(rownames(dag)))stop("dag is not a named matrix")
                        
  if(Reduce("|", names(data.dists)!=names(data.dists[names(data.df)])) | Reduce("|", names(data.dists)!=names(data.dists[colnames(dag)]))){
    stop("data.dists, data.df and dag do not have the same names or the same names' order")
  }
  
  if((!is.null(adj.vars) & !is.null(cor.vars)) & !(is.null(cor.vars[adj.vars]))){stop("cor.vars contains adj.vars, please remove them")}
  
  ## standardize gaussian variables to zero mean and sd=1
  if(centre && !is.null(data.dists=="gaussian")){## have at least one gaussian variable
    for(i in names(data.dists)[(data.dists=="gaussian")]){data.df[,i] <- (data.df[,i]-mean(data.df[,i]))/sd(data.df[,i]);}
  }
  
  ##formatting 
  for(i in 1:n){
    if(data.dists[[i]]=="binomial" & class(data.df[,i])!="numeric"){
      data.df[,i] <- as.numeric(factor(data.df[,i]))-1
    }
      if(data.dists[[i]]=="multinomial"){
        data.df[,i] <- factor(data.df[,i])
      }
  }
  
  ###----------------###
  ###start adjustment###
  ###----------------###
  
  if(!is.null(adj.vars)){
  
  if(is.null(cor.vars)){
    cor.vars <- colnames(data.df)
    cor.vars <- cor.vars[-adj.vars]}
  
  dag[cor.vars,adj.vars] <- 1
  }
  
    coef.out <- list()
    out <- list()
    ll.out.tmp <- list()
    aic.out.tmp <- list()
    bic.out.tmp <- list()
    mdl.out.tmp <- list()
    sse.out <- list()
    mse.out <- list()
    var.out <- list()
    pvalue <- list()
    
    #unpacking the multinomial variables in the data.df
    data.df.multi <- NULL
    
    for(i in 1:n){
      if(data.dists[[i]] %in% c("binomial", "poisson", "gaussian")){
        data.df.multi <- as.data.frame(cbind(data.df.multi,data.df[,i]))
        colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
      }else{
        tmp <- model.matrix(~-1+factor(data.df[,i]))
        colnames(tmp) <- paste0(names(data.df)[i], levels(factor(data.df[,i])))
        data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))}
    }
    
    #extend dag to multinomial variables
    repetition.multi <- vector(length = n)
    
    for(i in 1:n){
      if(data.dists[[i]] %in% c("binomial", "poisson", "gaussian")){
        repetition.multi[i] <- 1
      }else{
          repetition.multi[i] <- nlevels(data.df[,i])
          }
    }
    
    dag.multi <- dag[,rep(1:n, repetition.multi)]
    
    ##-----------------------------
    ##start loop for the regression
    ##-----------------------------
    
    if(verbose){cat("Start estimation loop.")}
    
    for(i in 1:length(dag[1,])){
      
      Y <- data.matrix(data.df[,i])
      
      if("multinomial" %in% data.dists[as.logical(dag[i,])]){
        X <- data.matrix(data.df.multi[,as.logical(dag.multi[i,])])
      }else{
        X <- data.matrix(cbind(rep(1,length(data.df[,1])),data.df[,as.logical(dag[i,])]))
      }
      
      num.na <- 0
      if(qr(X)$rank/ncol(X)!=1 & as.character(data.dists[i])=="binomial"){
        
        Y <- as.numeric(as.character(Y))
        
        fit <- tryCatch(irls_binomial_cpp_br(A = X, b = Y, maxit = maxit,tol = tol),error=function(e){
          while((qr(X)$rank/ncol(X))!=1){
            X <- X[,-1]
            num.na <- num.na+1
            if(is.null(ncol(X))) X <- as.matrix(X)
          }
          list(irls_binomial_cpp_br(A = X, b = Y, maxit = maxit,tol = tol),num.na)
        })
        num.na <- fit[[2]]
        fit <- fit[[1]]
      }else{

            switch(as.character(data.dists[i]),
             "binomial"={
               
              Y <- as.numeric(as.character(Y))
              
              #fit <- irls_binomial_cpp(A = X, b = Y, maxit = maxit,tol = tol)
              
              fit <- tryCatch(irls_binomial_cpp(A = X, b = Y, maxit = maxit,tol = tol),error=function(e){
                
                irls_binomial_cpp_br(A = X, b = Y, maxit = maxit,tol = tol)
              })
              
              if(is.na(sum(fit[[1]])))fit <- irls_binomial_cpp_br(A = X, b = Y, maxit = maxit,tol = tol)
              
             },
             "gaussian"={
               fit <- irls_gaussian_cpp(A = X, b = Y, maxit = maxit,tol = tol)
             },
             "poisson"={
               
               fit <- irls_poisson_cpp(A = X, b = Y, maxit = maxit,tol = tol)
             },
             "multinomial"={
               
               tmp <- multinom(formula = Y~-1+X,Hess = FALSE,trace=FALSE)
               
               #output
               fit <- list()
               fit$coefficients <- as.matrix(as.vector(coef(tmp)))
               fit$names.coef <- row.names((coef(tmp)))
               fit$loglik <- - tmp$value
               edf <- ifelse(length(tmp$lev) == 2L, 1, length(tmp$lev) - 1) * qr(X)$rank
               fit$aic <- 2 * tmp$value + 2 * edf
               fit$bic <- 2 * tmp$value + edf * log(nobs)
               fit$sse <- sum(residuals(tmp)^2)
               fit$var.out <- as.matrix(as.vector(summary(tmp)$standard.errors))
               
             }
             
      )}
      
      
      #fit$var.out <- t(fit$var.out)
      
      ll.out.tmp[[paste(names(data.dists[i]))]] <- fit$loglik
      aic.out.tmp[[paste(names(data.dists[i]))]] <- fit$aic
      bic.out.tmp[[paste(names(data.dists[i]))]] <- fit$bic
      mdl.out.tmp[[paste(names(data.dists[i]))]] <- fit$bic + (1 + sum(dag.multi[i,] - num.na)) * log(n)
      
      sse.out[[i]] <- fit$sse
      
      deg.freedom <- (length(data.df[,1])-(sum(dag.multi[i,])+1))
      
      #var.out[[i]] <- X
      mse.out[[i]] <- (fit$sse/deg.freedom)
      
      switch(as.character(data.dists[i]),
             "gaussian"={
               coef.out[[i]] <- matrix(data = c(rep(NA,num.na),fit$coefficients),nrow = 1)
               
               tmp.catch <- tryCatch(expr = solve(t(X) %*% X),error = function(e) NaN)
               if(is.nan(tmp.catch[1])){
                 var.out[[i]] <- matrix(data = rep(NaN,length.out=ncol(X)),nrow = 1)
                 pvalue[[i]] <- rep(1,length.out=ncol(X))
               }else{
                 var.out[[i]] <- matrix(data = sqrt(diag(mse.out[[i]]*tmp.catch)),nrow = 1)
                 var.out[[i]] <- matrix(data = c(rep(NA,num.na),var.out[[i]]),nrow = 1)
                 pvalue[[i]] <- 2*pt(-abs(coef.out[[i]]/var.out[[i]]),deg.freedom)
               }
               
             },
             "binomial"={
               coef.out[[i]] <- matrix(data = c(rep(NA,num.na),fit$coefficients),nrow = 1)
               #fit$coefficients <- rbind(as.matrix(rep(NA,num.na)),fit$coefficients)
               #coef.out[[i]] <- t(fit$coefficients)
               
               var.out[[i]] <- tryCatch(as.matrix(sqrt(diag(solve(fit$varcov)))),error=function(e){
                 
                 as.matrix(sqrt(svd(fit$varcov)$d))
               })
               
                 #var.out[[i]] <- as.matrix(sqrt(diag(solve(fit$varcov))))
                 var.out[[i]] <- matrix(data = c(rep(NA,num.na),var.out[[i]]),nrow = 1)
                 pvalue[[i]] <- 2*pnorm(-abs(coef.out[[i]]/var.out[[i]]))
             },
             "poisson"={
               fit$coefficients <- cbind(t(rep(NA,num.na)),t(fit$coefficients))
               coef.out[[i]] <- fit$coefficients
               
               var.out[[i]] <- as.matrix(sqrt(diag(solve(fit$varcov))))
               var.out[[i]] <- matrix(data = cbind(t(rep(NA,num.na)),t(var.out[[i]])),nrow = 1)
               pvalue[[i]] <- 2*pnorm(-abs(coef.out[[i]]/var.out[[i]]))
             },
             "multinomial"={
               fit$coefficients <- rbind(as.matrix(rep(NA,num.na*length(fit$names.coef))),fit$coefficients)
               coef.out[[i]] <- t(fit$coefficients)
               
               var.out[[i]] <- matrix(data = rbind(as.matrix(rep(NA,num.na*length(fit$names.coef))),fit$var.out),nrow = 1)
               pvalue[[i]] <- 2*pnorm(-abs(coef.out[[i]]/var.out[[i]]))
             }
      )
      
      #naming
      #var.out[[i]] <- t(matrix(var.out[[i]]))
      pvalue[[i]] <- matrix(data = pvalue[[i]],nrow = 1)
      
      ##first if response is multinomial
      
      if(data.dists[[i]]!="multinomial"){

        if("multinomial" %in% data.dists[as.logical(dag[i,])]){
          colnames(coef.out[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
          colnames(var.out[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
          colnames(pvalue[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
        }else{
          colnames(coef.out[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
          colnames(var.out[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
          colnames(pvalue[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
        }
      }
        if(data.dists[[i]]=="multinomial"){
          if("multinomial" %in% data.dists[as.logical(dag[i,])]){
            colnames(coef.out[[i]]) <- c(as.vector(outer(names(data.df.multi)[as.logical(dag.multi[i,])], fit$names.coef, paste, sep=".")))
            colnames(var.out[[i]]) <- c(as.vector(outer(names(data.df.multi)[as.logical(dag.multi[i,])], fit$names.coef, paste, sep=".")))
            colnames(pvalue[[i]]) <- c(as.vector(outer(names(data.df.multi)[as.logical(dag.multi[i,])], fit$names.coef, paste, sep=".")))
          }else{
            colnames(coef.out[[i]]) <- c(paste(names(data.df)[i],"|intercept.",fit$names.coef,sep = ""),as.vector(outer(names(data.df)[as.logical(dag[i,])], fit$names.coef, paste, sep=".")))
            colnames(var.out[[i]]) <- c(paste(names(data.df)[i],"|intercept.",fit$names.coef,sep = ""),as.vector(outer(names(data.df)[as.logical(dag[i,])], fit$names.coef, paste, sep=".")))
            colnames(pvalue[[i]]) <- c(paste(names(data.df)[i],"|intercept.",fit$names.coef,sep = ""),as.vector(outer(names(data.df)[as.logical(dag[i,])], fit$names.coef, paste, sep=".")))
          }
        }
      
      
      
      if((("multinomial" %in% data.dists[as.logical(dag[i,])]) & (data.dists[[i]]!="multinomial"))){
        colnames(coef.out[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
        colnames(var.out[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
        colnames(pvalue[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
      }
      #ok
      if(!("multinomial" %in% data.dists)){
        colnames(coef.out[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
        colnames(var.out[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
        colnames(pvalue[[i]]) <- c(paste(names(data.df)[i],"|intercept",sep = ""),names(data.df)[as.logical(dag[i,])])
      }
      
      # if(("multinomial" %in% data.dists[as.logical(dag[i,])]) & (data.dists[[i]]=="multinomial")){
      #   colnames(coef.out[[i]]) <- c(paste(names(data.df.multi)[as.logical(dag.multi[i,])],fit$names.coef))
      #   colnames(var.out[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
      #   colnames(pvalue[[i]]) <- c(names(data.df.multi)[as.logical(dag.multi[i,])])
      # }
      # if(data.dists[[i]]=="multinomial" & !("multinomial" %in% data.dists[as.logical(dag[i,])])){
      # 
      # }
    }#EOF loop Regression

    names(coef.out) <- names(data.dists)
    names(var.out) <- names(data.dists)
    names(pvalue) <- names(data.dists)
    
    #mlik
    #lapply(l, function(x) replace(x, which(is.na(x)), 0))
    out[["method"]] <- "mle"
    
    out[["abnDag"]] <- createAbnDag(dag, data.df = data.df, data.dists = data.dists)
     
    out[["mliknode"]] <- ll.out.tmp
    out[["mlik"]] <- Reduce("+",ll.out.tmp)
    out[["aicnode"]] <- aic.out.tmp
    out[["aic"]] <- Reduce("+", aic.out.tmp)
    out[["bicnode"]] <- bic.out.tmp
    out[["bic"]] <- Reduce("+", bic.out.tmp)
    out[["mdlnode"]] <- mdl.out.tmp
    out[["mdl"]] <- Reduce("+", mdl.out.tmp)
    out[["coef"]] <- (coef.out)
  
    #p value
    out[["Stderror"]] <- var.out
    out[["pvalue"]] <- pvalue
    #out[["mse"]] <- mse.out
    
    class(out) <- c("abnFit")
    
    return(out)
  }
