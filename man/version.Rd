\name{version}
\alias{version}
\alias{abn.version}
\alias{abn.Version}
\title{abn Version Information}
\description{
  \code{abn.version} is a variable (class \code{simple.list}) holding
  detailed information about the version of \code{abn} loaded.

  \code{abn.Version()} provides detailed information about the running version of \pkg{abn}
     or the \pkg{abn} components.
}
\usage{
abn.Version(what=c("abn","system"))
abn.version
}
\arguments{
  \item{what}{detailed information about the version of \pkg{abn}
     or a summary of the \pkg{abn} components.}
     }
\value{
  \code{abn.Version()} is a list with character-string components
  \item{R}{\code{R.version.string}}
  \item{abn}{essentially \code{abn.version$version.string}}
  \item{GSL, JAGS, INLA}{version numbers thereof}

  \code{abn.version} is a list with character-string components
  \item{status}{the status of the version (e.g., \code{"beta"})}
  \item{major}{the major version number}
  \item{minor}{the minor version number}
  \item{year}{the year the version was released}
  \item{month}{the month the version was released}
  \item{day}{the day the version was released}
  \item{version.string}{a \code{character} string concatenating
     the info above, useful for plotting, etc.}

  \code{abn.version} is a list of class \code{"simple.list"} which
     has a \code{print} method.
}
% \references{}
\seealso{See the R counterparts  \code{\link[base]{R.version}}.
}
\author{Reinhard Furrer}
\examples{
abn.version$version.string

\dontrun{
abn.Version("system")
}
}
\keyword{environment}
\keyword{sysdata}
\keyword{programming}
