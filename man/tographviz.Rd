%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{tographviz}
\alias{tographviz}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Convert a dag into graphviz format}

\description{Given a matrix defining a DAG create a text file suitable for plotting with graphviz}


\usage{tographviz(dag.m,data.df,data.dists,group.var,outfile)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag.m}{a matrix defining a DAG}
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{data.dists}{a list with named arguments matching the names of the data frame which gives the distribution family for each variable. See \code{fitabn} for details.} 
  \item{group.var}{only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership). See \code{fitabn} for details.}
  \item{outfile}{a character string giving the filename which will contain the graphviz graph}
}

\details{Graphviz - \url{www.graphviz.org} is visualisation software developed by AT&T and freely available. This function creates a text representation of the DAG so this can be plotted using graphviz. Graphviz is available as an R package, Rgraphviz, through the Bioconductor project \url{http://www.bioconductor.org/} (and requires a working installation of graphviz). Binary nodes will appear as squares, Gaussian as ovals and Poisson nodes as diamonds in the resulting graphviz network diagram. There are many other shapes possible for nodes and numerous other visual enhancements - see online graphviz documentation. Bespoke refinements can be added by editing the raw outfile produced. For full manual editing, particularly of the layout, or adding annotations, one easy solution is to convert a postscript format graph (produced in graphviz using the -Tps switch) into a vector format using a tool such as pstoedit \url{www.pstoedit.net}, and then edit using a vector drawing tool like xfig. This can then be resaved as postscript or pdf thus retaining full vector quality} 


\value{nothing is returned}
\references{

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.r-bayesian-networks.org}}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\examples{
\dontrun{

mydat<-ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

## setup distribution list for each node
mydists<-list(b1="binomial",
              b2="binomial",
              b3="binomial",
              g1="gaussian",
              b4="binomial",
              p2="poisson",
              p4="poisson"
             );
## specify DAG model 
mydag<-matrix(data=c(
                     0,1,0,0,1,0,0, # 
                     0,0,0,0,0,0,0, #
                     0,1,0,0,1,0,0, # 
                     1,0,0,0,0,0,1, # 
                     0,0,0,0,0,0,0, #
                     0,0,0,1,0,0,0, #
                     0,0,0,0,1,0,0  #
                     ), byrow=TRUE,ncol=7);
colnames(mydag)<-rownames(mydag)<-names(mydat);

## create file for processing with graphviz
tographviz(dag.m=mydag,data.df=mydat,data.dists=mydists,outfile="graph.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graph.pdf graph.dot")
system("evince graph.pdf");

## example using data with a group variable.

## model where b1<-b2
mydag<-matrix(data=c(
                     0,1, # b1
                     0,0  # b2
                     ), byrow=TRUE,ncol=2);

colnames(mydag)<-rownames(mydag)<-names(ex3.dag.data[,c(1,2)]);
## specific distributions
mydists<-list(b1="binomial",
              b2="binomial"
             );

## create file for processing with graphviz
tographviz(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
           group.var="group",outfile="graph.dot");
## and then process using graphviz tools e.g. on linux
system("dot -Tpdf -o graph.pdf graph.dot")
system("evince graph.pdf");
}

}

\keyword{graphviz}
