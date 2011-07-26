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


\usage{tographviz(dag,outfile)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag}{a matrix defining a DAG}
  \item{outfile}{a character string giving the filename which will contain the graphviz graph}
}

\details{Graphviz - \url{www.graphviz.org} is visualisation software developed by AT&T and freely available. This function creates a text representation of the DAG so this can be plotted using graphviz. Graphviz is available as an R package, Rgraphviz, through the Bioconductor project \url{http://www.bioconductor.org/} }  



\value{nothing is returned}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\examples{
dag<-matrix(c(0,0,0,0,0,0,0,0,0,0, ## var1
                 0,0,0,0,1,0,0,0,0,0, ## var2<-var5
                 0,0,0,0,1,0,0,0,0,0, ## var3<-var5
                 0,0,0,0,0,0,0,0,0,0, ## var4
                 0,0,0,0,0,0,0,0,0,0, ## var5
                 1,0,0,0,0,0,0,0,0,0, ## var6<-var1
                 0,1,0,0,0,0,0,0,0,0, ## var7<-var2
                 1,0,1,0,0,0,0,0,0,0, ## var8<-var1:var3
                 0,0,0,0,1,0,0,0,0,0, ## var9<-var5
                 0,0,0,1,0,0,0,0,0,0),## var10<-var4
                      byrow=TRUE,ncol=10);

colnames(dag)<-c("var1","var2","var3","var4","var5",
                 "var6","var7","var8","var9","var10");
rownames(dag)<-colnames(dag);
\dontrun{tographviz(dag,outfile="mygraph.dot");}
\dontrun{system("dot -Tps -o mygraph.ps mygraph.dot");}
## dot is part of graphviz on unix

}

\keyword{graphviz}
