\name{sbd}
\alias{sbd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates sample bivariate data.
}
\description{
  This function generates a sample bivariate data set. 
}
\usage{
sbd(func,min,max,n,Rsq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
    a user supplied function of one variable, y = func(x), near which data is generated.
  }
  \item{min}{
    min value for the domain of func
  }
  \item{max}{
    max value for the domain of func
  }
  \item{n}{
    number of sample points to generate
  }
  \item{Rsq}{
    coefficient of determination for the data set
  }
}
\details{
  If func is NULL then a normal bivariate data set of n samples is generated
  with correlation coefficient sqrt(Rsq).
  If func is passed by the user then n sample points are scattered about y=func(x)
  with variance governed by the Rsq parameter
}
\value{
  Returns an n x 2 bivariate data set
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
  See examples below on how to set up user defined functions
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
}

\examples{
    f <- function(x,name="Sinusoidal",def="y = 1 + sin(x)"){
      return(1 + sin(x))
    }
    d <- sbd(f,min=-2*pi,max=2*pi,n=500,Rsq=0.9)
    ma(d)$A
    plot(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{coefficient of determination}
