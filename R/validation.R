#
# Ported function to check the various classes. Mostly refers to abn-internal.R


create_abnDag <- function( dag, data.df=NULL, data.dists=NULL, ...) {
  
  if (!is.null(data.dists))   
    data.dists <- validate_dists(data.dists, returnDists=TRUE)
  
  if (is.null( data.df)) {
    dag <- validate_abnDag(dag, data.df=data.dists, returnDag=TRUE)
  }else{
    dag <- validate_abnDag(dag, data.df=data.df, returnDag=TRUE)
  }
  
  out <- list( dag=dag, data.df=data.df, data.dists=data.dists) 
  class( out) <- "abnDag"
  
  return( out)
  
}

validate_dists <- function(data.dists, returnDists=TRUE,...) {
  
  name <- names(data.dists)
  if (is.null(name)) stop("Node distribution has to be a named object.")
  if( is.list( data.dists))       data.dists <- unlist( data.dists)
  
  choices <- c("poisson","binomial","gaussian","multinomial")
  data.dists <- choices[pmatch(tolower(data.dists ), choices, duplicates.ok=TRUE)]
  if (any(is.na(data.dists ))) stop("Incorrectly specified node distribution.")
  names(data.dists ) <- name
  
  
  if( returnDists) return( as.list( data.dists)) else return( TRUE)
  
}

validate_abnDag <- function( dag, data.df=NULL, returnDag=TRUE, ...) {
  
  # dag is either a formula or a matrix 
  
  # we already have a valid container. can beused to extract...  
  if (inherits(x = dag, what = "abnDag"))  dag <- dag$dag
  
  
  # case of formula
  if  (inherits(x = dag, what = "formula")) {
    if (is.null( data.df))
      stop( 'DAG specification with formula requires a named data frame or named vector')
    
    name <- if ( is.matrix( data.df)) colnames( data.df) else names( data.df)
    if (is.null( name)) 
      stop( 'Improperly named object "data.df" for DAG specification')
    
    dag <- formula.abn(f = dag, name = name)
  }   # proceed checking!!
  
  # case of matrix
  
  if ( is.matrix( dag)) {
    dimm <- dim( dag) 
    if (dimm[1] != dimm[2])   stop("DAG matrix is not square")
    if (any(diag(dag)!=0))  stop("DAG matrix contains trivial cycles (nonzero values on diagonal)")
    
    if (!is.null(data.df))  {    # if data.df given we take over the names.
      name <- if ( is.matrix( data.df)) colnames( data.df) else names( data.df)
      
      if(length(name) != dimm[1])  stop("DAG matrix not coherent with names")
      colnames(dag) <- rownames(dag) <- name
    } else {
      if (any(colnames(dag)!=rownames(dag)))  stop("DAG matrix with incoherent row-/colnames")
    }
    
    res <- .Call("checkforcycles", as.integer(dag), dimm[1], PACKAGE = "abn")
    
    if( returnDag) return( dag) else return( TRUE)
  }   else {
    stop( "DAG specification with should be via formula or matrix")
  }
}