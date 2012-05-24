\name{spa}
\alias{spa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Semi-Partial Association (for a set of 3 variables) 
}
\description{
Computes the semi-partial association between a response variable and 
an explanatory variable, after controlling for a control variable. 
}
\usage{
spa(Y,X,C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  the response variable, a vector or column from a dataset
  }
  \item{X}{
  the explanatory variable, a vector or column from a dataset
  }
  \item{C}{
  the control variable, a vector or column from a dataset
  }
}

\details{
  A semi-partial association (possibly nonlinear) is computed via:
  
  \code{ ma(cbind(C,X,Y))$A - ma(cbind(C,Y))$A }.
  
  Inspired by the linear semi-partial correlation given by:
  
  \code{ spcor.test(Y,X,C) }
  
  from the \code{ppcor} package.
  
}
  

\value{
  Returns a real number in the range [0,1].
}

\references{
coming soon
}

\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}

\note{
The parameters Y, X and C must be vectors of the same length.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
}

\examples{
    # Example showing similar behavior of linear and non-linear 
    # semi-partial associations when the data is truly linear
    # Y is the response variable, X is the explanatory variable, and C is the control
    # k is a factor that controls how much of Y is taken from X and how much from C
    k <- 0.5
    C = rnorm(1000)
    X = rnorm(1000)
    Y = C*(1-k) + X*k + rnorm(1000)*0.1

    # if you have rgl you can plot the data
    # plot3d(X,C,Y)
    
    # compute the semi-partial association
    spa(Y,X,C)
    
    # if you have ppcor then you can compute the linear semi-partial correlation as well
    # spcor.test(Y,X,C)$estimate^2

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-partial }
\keyword{ correlation }
\keyword{ association }
