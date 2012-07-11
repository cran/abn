%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{getmarginal}
\alias{getmarginal}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Estimation of marginal posterior density in an additive Bayesian network model}

\description{Estimates the marginal posterior density for a single regression coefficient in an additive Bayesian network model}

\usage{
getmarginal(data.df,dag.m,whichnode=NULL, whichvar="constant", 
                        hyper.params=list("mean"=c(0),
                                          "sd"=c(sqrt(1000)),
                                          "shape"=c(0.001),
                                          "scale"=c(1/0.001)),
            verbose=FALSE,post.x=seq(-0.1,0.1,len=100),max.iters=100,
            epsabs=1e-7,std=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network, will be coerced to factors if necessary}
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph, see details for format. Note that colnames and rownames must be set}
  \item{whichnode}{the name of the node - response variable}
  \item{whichvar}{the name of the explanatory variable in whichnode whose posterior density is required. For the intercept term use ``constant'', for the precision parameter for a Gaussian node use ``precision''}
  \item{hyper.params}{prior means and standard deviations in Gaussian prior distribution for regression coefficients, see details}
  \item{verbose}{logical, if TRUE then extra information is printed to the terminal} 
  \item{post.x}{the set of points in the domain at which to evaluate the posterior density}
  \item{max.iters}{maximum number of iterations used to find root in Laplace approximation}
  \item{epsabs}{minimum absolute error of root estimate in Laplace approximation}
\item{std}{logical, if TRUE then each variable which is not a factor is standardised to zero mean and standard deviation one} 
}

\details{
  The procedure \code{getmarginal} estimates the posterior density for a regression coefficient within an additive Bayesian network model and is similar to \code{fitabn} but extracts the posterior density rather than the log marginal likelihood. The user must supply sensible values for post.x, e.g. a range which is sufficient to include both left and right tail of the density. As fixed grid is used and so using excessively wide ranges will be inefficient or result in very poor accuracy. Sensible values for post.x can easily be found by some trial and error, and will generally be close to the origin due to the use of a logit link or the standardization of Gaussian variables. 
  
  In the network structure definition, dag.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.
  

}

\value{A matrix with two columns called x and f, x is the value of the variable and f is the corresponding density value.}


\references{Tierney, L. and Kadane, J. B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities.
Journal of the American Statistical Association, Vol. 81, No. 393, pp. 82-86.


  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitabn}}}

\examples{
#specific a DAG model
mydag<-matrix(c(
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#1
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#2  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#3  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#4  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#5  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#6  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#7  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#8  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#9  
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#10 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#11 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#12 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#13 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#14 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#15 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#16 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#17 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#18 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#19 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#20 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#21 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#22 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#23 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#24 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#25 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#26 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#27 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#28 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#29 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#30
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#31
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,#32 
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 #33
                                               ),byrow=TRUE,ncol=33); 
colnames(mydag)<-rownames(mydag)<-names(var33);#set names
## now fit the model defined in mydag - full independence
# define a model with many independencies
mydag[2,1]<-1;
mydag[4,3]<-1;
mydag[6,4]<-1; mydag[6,7]<-1;
mydag[5,6]<-1;
mydag[7,8]<-1;  
mydag[8,9]<-1;
mydag[9,10]<-1;
mydag[11,10]<-1; mydag[11,12]<-1; mydag[11,19]<-1;
mydag[14,13]<-1;
mydag[17,16]<-1;mydag[17,20]<-1;
mydag[15,14]<-1; mydag[15,21]<-1;
mydag[18,20]<-1;
mydag[19,20]<-1;
mydag[21,20]<-1;
mydag[22,21]<-1;
mydag[23,21]<-1;
mydag[24,23]<-1;
mydag[25,23]<-1; mydag[25,26]<-1;
mydag[26,20]<-1;
mydag[33,31]<-1;
mydag[33,31]<-1;
mydag[32,21]<-1; mydag[32,31]<-1;mydag[32,29]<-1;    
mydag[30,29]<-1;
mydag[28,27]<-1; mydag[28,29]<-1;mydag[28,31]<-1;
marg1<-getmarginal(data.df=var33,
                   dag.m=mydag,
                   whichnode="v33",
                   whichvar="constant",#this is the intercept
                    hyper.params=list(
                     mean=rep(0,dim(var33)[2]+1),
                     sd=rep(sqrt(1000),dim(var33)[2]+1), 
                     shape=rep(0.001,16),## 16 Gaussian nodes
                     scale=rep(1/0.001,16)## 16 Gaussian nodes
                     ),
                     post.x=seq(from=-1.5,to=-0.5,len=1000),
                     verbose=TRUE);
cum.marg<-cumsum(marg1[,"f"])/sum(marg1[,"f"])
marg<-cbind(marg1,cum.marg);
marg[which(marg[,3]>0.5)[1]];#approx.median

}

\keyword{models}
