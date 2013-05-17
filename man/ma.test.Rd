\name{ma.test}
\alias{ma.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
p-value for an association score.
}
\description{
  This function computes the p-value for an association score using 
  interpolation from a precomputed table 
  (fast - only available for the bivariate case) or using a permutation test.
}
\usage{
ma.test(d, maStruct, permPartition, lookup, reps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
    the data set
  }
  \item{maStruct}{
    the structure returned by a call to \code{\link{ma}}.
    Must be provided. No default.
  }
  \item{permPartition}{
    a permutation partition for the data set depending on 
    the statistical question being asked.
    Defaults to the partition returned by the call to \code{\link{ma}}.
  }
  \item{lookup}{
    logical, if TRUE then a p-value is interpolated from a table of 
    precomputed p-values. 
    If FALSE then a monte carlo simulation is carried out. 
    Table lookup is only accurate for associations on bivariate data 
    with number of samples <= 400. By default, lookup is TRUE.
  }
  \item{reps}{
    no of repeats for the monte carlo simulation.
    defaults to 1000.
  }
}
\details{
  Computes a p value for an association from the likelihood ratio statistic
  computed by \code{\link{ma}}. Uses interpolation by default for bivariate associations,
  and permutation testing for multivariate associations.
}
\value{
  Returns a real number in the range [0,1].
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
}

\examples{
    d <- shpd(n=100,m=2,Rsq=0.2)
    aScore <- ma(d)
    ma.test(d, aScore)
    # or set lookup to FALSE if you have computing time available
    # ma.test(d, aScore, lookup=FALSE, reps=10000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ p value }
