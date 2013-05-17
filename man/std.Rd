\name{std}
\alias{std}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates sample trivariate data set.
}
\description{
  This function generates a sample trivariate data set. 
}
\usage{
std(func,xMin,xMax,yMin,yMax,n,Rsq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
    a user supplied function of two variables, z = func(x,y), near which data is generated
  }
  \item{xMin}{
    min value for the x domain of func
  }
  \item{xMax}{
    max value for the x domain of func
  }
  \item{yMin}{
    min value for the y domain of func
  }
  \item{yMax}{
    max value for the y domain of func
  }
  \item{n}{
    number of sample points to generate
  }
  \item{Rsq}{
    coefficient of determination for the data set
  }
}
\details{
  If func is NULL then a normal trivariate data set of n samples is generated
  with correlation coefficients all set to sqrt(Rsq).
  If func is passed by the user then n sample points are scattered about z=func(x,y)
  with variance governed by the Rsq parameter.
}
\value{
  Returns an n x 3 trivariate data set.
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
  See examples below on how to set up user defined functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
\code{\link{sbd}}
}

\examples{
    f <- function(x,y,name="MexicanHat", def="z=(1-t^2)exp(-t^2div2), t^2=x^2+y^2"){
      t <- sqrt(x^2 + y^2)
      z <- (1.0 - t^2) * exp(- t * t / 2)
      return(z)
    }
    d <- std(f, xMin=-2, xMax=2, yMin=-2, yMax=2, n=500, Rsq=0.85)
    ma(d)$A
    # if you have rgl you can view the data set in 3D 
    # library("rgl")  
    # plot3d(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{coefficient of determination}