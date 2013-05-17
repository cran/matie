\name{fbvs}
\alias{fbvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forwards, backwards variable selection, 
picking variables to maximize explained variance.
}
\description{
Attempts to find the best explanatory set of variables to explain a single variable in
a data set. Iterates between adding the next best variable to the set and removing the variable (if any) whose exclusion maximizes the overall score.
}
\usage{
fbvs(dataSet,one,maxv,linear)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
  the n x m data frame representing n observations of m variables.
  }
  \item{one}{
  a string specifying the name of one variable in the dataset,
  for which the best explanatory set is required. Defaults to
  the name of the last variable in the dataset.
  }
  \item{maxv}{
  an integer limiting the maximum number of variables in the explanatory set.
  Defaults to m-1.
  }
  \item{linear}{
  a boolean flag which causes fbvs to use a linear model to estimate R^2
  instead of matie to estimate A when running the selection algorithm.
  Defaults to \code{FALSE}
  }
}
\details{
  Variable names are only added to the explanatory set if their inclusion
  results in an increase in the association measure. 
}
\value{
  Returns a list containing the following items:
  \item{one}{
  the name of the one variable that requires the explanatory set
  }
  \item{best}{
  the best set of explanatory variables
  }
  \item{Rsq}{
  an estimate for R^2 provided by the best set of explanatory variables
  }
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
The data set can be of any dimension
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
\code{\link{agram}}
}

\examples{
    # measure association for all pairs in a subrange of the baseball dataset
    data(baseballData)
    fbvs(baseballData,one="Salary")
    fbvs(baseballData,one="Salary",linear=TRUE)

    fbvs(baseballData,one="Salary",maxv=2)
    fbvs(baseballData,one="Salary",maxv=2,linear=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ association }
