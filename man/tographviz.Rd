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


\usage{tographviz(dag,data.df,outfile)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dag}{a matrix defining a DAG}
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{outfile}{a character string giving the filename which will contain the graphviz graph}
}

\details{Graphviz - \url{www.graphviz.org} is visualisation software developed by AT&T and freely available. This function creates a text representation of the DAG so this can be plotted using graphviz. Graphviz is available as an R package, Rgraphviz, through the Bioconductor project \url{http://www.bioconductor.org/} (and requires a working installation of graphviz). The data.df argument is used to determine which nodes are binary and which are gaussian, where the former will appear as squares in the resulting graphviz network diagram and the latter as ovals. There are many other shapes possible for nodes  and numerous other visual enhancements - see graphviz documentation - which can be added by editing the raw outfile produced. For full manual editing, particularly of the layout, one easy solution is to convert a postscript graph into a vector format using pstoedit \url{www.pstoedit.net} and edit using a vector drawing tool like xfig.} 


\value{nothing is returned}


\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\examples{

## see abn vignette for examples

}

\keyword{graphviz}
