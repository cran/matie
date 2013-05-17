\name{pd}
\alias{pd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot density (after an \code{ma} calculation).
}
\description{
  Optimal kernel widths output by \code{\link{ma}} are employed to 
  recompute the weighted joint distribution for two variables in
  a data set, and a contour plot for this distribution is drawn.
}
\usage{
pd(d,iv=1,jv=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
  an n x m data frame with m > 1.
  }
  \item{iv}{
  the column index of the independent variable
  }
    \item{jv}{
  the column index of the dependent variable
  }
}
\details{
  A data set of two variables is extracted from the user's data set and a 
  full distribution is calculated
  using weighted marginal and joint likelihoods. The optimal kernel sizes and 
  weighting are first computed via a call to \code{\link{ma}}.
}
\value{
  An n x n  distribution of weighted likelihoods is returned.
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
  Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
  The data set must contain at least 2 columns.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
}
\examples{
    f <- function(x,name="Sinusoidal",def="y = 1 + sin(x)"){
      return(1 + sin(x))
    }
    d <- sbd(f,min=-2*pi,max=2*pi,n=1000,Rsq=0.9)
    ma(d)$A
    dist <- pd(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ density }
\keyword{ contour }
% __ONLY ONE__ keyword per line
