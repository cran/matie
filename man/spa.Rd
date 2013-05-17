\name{spa}
\alias{spa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Semi-partial association (computes association while controlling for variables)
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
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
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
    # When the association between Y and X is entirely explained by C
    C <- rnorm(1000)
    X <- C+rnorm(1000)*0.5
    Y <- C+rnorm(1000)*0.5
    # See the relationship without controlling for C
    ma(cbind(X,Y))$A
    # See the relationship with C as a covariate (should be close to 0)
    spa(Y,X,C)
    # if you have ppcor then you can verify that 
    # the linear semi-partial correlation is similar, 
    # as these associations are all linear
    # spcor.test(Y,X,C)$estimate^2
    #
    # When the association between Y and X is only partially explained by C
    C <- rnorm(1000)
    X <- C+rnorm(1000)*0.5
    Y <- X+rnorm(1000)*0.5
    # See the relationship without controlling for C
    ma(cbind(X,Y))$A
    # See the relationship with C as a covariate 
    # (should be lower than the uncontrolled one, but not as low as 0)
    spa(Y,X,C)
    # if you have ppcor then you can verify that 
    # the linear semi-partial correlation is similar, 
    # as these associations are all linear
    # spcor.test(Y,X,C)$estimate^2
    #
    #   
    # if you have rgl you can plot the data
    # library(rgl)
    # plot3d(X,C,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ semi-partial }
\keyword{ correlation }
\keyword{ association }
