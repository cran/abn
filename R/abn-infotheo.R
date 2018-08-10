###############################################################################
## abn-infotheo.R ---
## Author          : Gilles Kratzer
## Document created   : 14/02/2016
## Last modification  : 14/02/2016
## Last modification  : 08/09/2017 Correction MI/Shanon entropy implementation
## Last modification  : 13/07/2018 cleaning + drag from varrank
###############################################################################


##---------------------------------------------------------------------------------------------------------
## Discretization of an arbitrary large number of variables joint variables depending on their distribution
## Implemented discretization methods:
## -user defined
## -fd
## -doane
## -sqrt
## -sturges
## -rice
## -scott
## -terrell-scott
## -cencov
##---------------------------------------------------------------------------------------------------------

discretization <- function(data.df = NULL, 
                           data.dists = NULL, 
                           discretization.method = "sturges", 
                           nb.states = FALSE){
  
  #tests:
  
  if(is.character(discretization.method)) discretization.method <- tolower(discretization.method)
  
  if(is.character(discretization.method)) discretization.method <- c("fd","doane","cencov","sturges","rice","scott","sqrt","terrell-scott")[pmatch(discretization.method,c("fd","doane","cencov","sturges","rice","scott","sqrt","terrell-scott"))]
  
  if(!(discretization.method %in% c("fd","doane","cencov","sturges","rice","scott","sqrt","terrell-scott") | is.numeric(discretization.method))){stop("Wrong definition of the discretization.method's parameter")}

  df.tab<-list()
  nb.states.l<-list()
  data.df<-as.data.frame(data.df)

  ##===============================================================
  ## Discrete cutoff provided
  ##===============================================================

  if(is.numeric(discretization.method)){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = discretization.method + 1)
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Freedman-Diaconis rule
  ##===============================================================

  if(discretization.method=="fd"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = (abs(range(data.df[,i])[2] - range(data.df[,i])[1]) /diff(range(data.df[,i]))/(2 * IQR(data.df[,i]) / length(data.df[,i])^(1/3)))+1)
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Doane s formula
  ##===============================================================

  if(discretization.method=="doane"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = 1 + log2(length(data.df[,i])) + log2( 1 + (abs(skewness(data.df[,i]))/sqrt((6*(length(data.df[,i])-2))/((length(data.df[,i])+1)*(length(data.df[,i])+3))))))
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Square root choice
  ##===============================================================

  if(discretization.method=="sqrt"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = ceiling(abs(sqrt(length(data.df[,i]))+1)))
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Sturges
  ##===============================================================

  if(discretization.method=="sturges"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = log(x = length(data.df[,i]),base = 2)+2)
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Rice Rule
  ##===============================================================

  if(discretization.method=="rice"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = 2*length(data.df[,i])^(1/3)+1)
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}

  ##===============================================================
  ##Scott rule
  ##===============================================================

  if(discretization.method=="scott"){

    for(i in 1:length(data.dists)){

      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = (abs(range(data.df[,i])[2] - range(data.df[,i])[1]) /diff(range(data.df[,i]))/(3.5 * sd(data.df[,i]) / length(data.df[,i])^(1/3)))+1)
      }

      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}
  ##===============================================================
  ##Terrell-Scott
  ##===============================================================
  
  if(discretization.method=="terrell-scott"){
    
    for(i in 1:length(data.dists)){
      
      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = ceiling(abs(2*nobs)^(1/3)))
      }
      
      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}
  
  
  
  ##===============================================================
  ##cencov
  ##===============================================================
  
  if(discretization.method=="cencov"){
    
    for(i in 1:length(data.dists)){
      
      if(data.dists[i]=="binomial" | data.dists[i]=="multinomial"){
        cut.index1 <- seq(from = range(as.numeric(data.df[,i]))[1], to = range(as.numeric(data.df[,i]))[2], length.out = length(levels(data.df[,i]))+1)
      }
      else{
        cut.index1 <- seq(from = range(data.df[,i])[1], to = range(data.df[,i])[2], length.out = ceiling(abs(nobs^(1/3))))
      }
      
      df.tab[[i]] <- as.factor(cut(as.numeric(data.df[,i]), breaks = cut.index1, include.lowest = TRUE))
      nb.states.l[[i]] <- length(cut.index1)
    }}
  
  



if(nb.states){
  return(list("df"=table(df.tab),"nb.states"=nb.states.l))
}else{
    return(table(df.tab))
  }
}

##-------------------------------------------------------------------------
## Mutual Information
## Function that returns MI or MI in percentage
##-------------------------------------------------------------------------

miData <- function(freqs.table, method = c("mi.raw", "mi.raw.pc")){

  #normalization
  freqs.joint <- prop.table(freqs.table)

  #xyz ...
  for(i in 1:length(dim(freqs.joint))){
    if(i==1){freqs.xy<-margin.table(x = freqs.table,margin = 1)
    }else{
      freqs.xy <- outer(X = freqs.xy,Y = margin.table(x = freqs.table,margin = i),FUN = "*")}
  }

  freqs.xy <- freqs.xy/sum(freqs.xy)

  #computing log part
  log.part <- ifelse(freqs.joint > 0, log(freqs.joint/freqs.xy), 0)

  #computing MI
  mi <- sum(freqs.joint * log.part)

  if(method == "mi.raw"){return(mi)}

  if(method == "mi.raw.pc"){
    tmp<-ifelse((1/margin.table(x = freqs.table,margin = 1))<Inf,log(1/margin.table(x = freqs.table,margin = 1),base = 2),0)
    return(-(mi*100)/sum(margin.table(x = freqs.table,margin = 1)*tmp))
    #
  }

  if(FALSE){

  #normalization
  freqs.joint <- as.matrix(freqs.table/sum(freqs.table))

  #xy
  freqs.x <- rowSums(freqs.table)
  freqs.y <- colSums(freqs.table)
  freqs.xy <- freqs.x %o% freqs.y
  freqs.xy <- freqs.xy/sum(freqs.xy)

  #computing log part
  log.part <- ifelse(freqs.joint > 0, log(freqs.joint/freqs.xy), 0)

  #computing MI
  mi <- sum(freqs.joint * log.part)

  if(method == "mi.raw"){return(mi)}

  if(method == "mi.raw.pc"){
    tmp<-ifelse((1/freqs.x)<Inf,log(1/freqs.x,base = 2),0)
    return(-(mi*100)/sum(freqs.x*tmp))
    #
  }
}

}

##-------------------------------------------------------------------------
## Shanon Entropy
##-------------------------------------------------------------------------

entropyData <- function(freqs.table){

  #normalization
  freqs.joint <- prop.table(freqs.table)

  #computing log part
  log.part <- ifelse(freqs.joint > 0, log(freqs.joint,base = 2), 0)

  #computing entropy
  entropy.out <- sum(freqs.joint * log.part)

  return(-entropy.out)

}

##EOF
