\name{tap}
\alias{tap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test all pairs (computes association score for all pairs of variables).
}
\description{
Tries all pairs in an n x m user supplied data set (n observations of m variables).
Returns a square matrix of association measures. 
}
\usage{
tap(dataSet,one)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
  the n x m data frame representing n observations of m variables.
  }
  \item{one}{
  a string specifying the name of one variable in the dataset,
  for testing all pairs against one.
  }
}
\details{
  All pairs of variables in a user supplied data set are measured for A association. 
  Returns an \code{m x m} data frame ( or \code{m-1 x m-1} when one is set) 
  storing A values for each pair of variables in the all pairs test. 
}
\value{
  Returns a square data frame.
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
The data set can be of any dimension.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
}

\examples{
    # measure association for all pairs in a subrange of the baseball dataset
    data(baseballData)
    ap <- tap(baseballData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ association }
