\name{ma.nl}
\alias{ma.nl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculates nonlinear preportion of variance between one variable
  and a group of variables.
}
\description{
  This function estimates nonlinear proportion of variance between one variable
  and a group of variables after linear association between the variable
  and the group has been removed.
}
\usage{
ma.nl(Y, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
    A vector or a one column data frame.
  }
  \item{X}{
    a group of vectors or a data frame with the same number of samples as in Y
  }
}
\details{
  A linear model, \code{Y ~ X}, is constructed and \code{ma}
  is used to compute R^2 between Y and X.  
}
\value{
  Returns a list of real numbers:
  \item{Rsq}{
    linear association, the value of R^2 due to the linear model \code{Y ~ X}.
  }
  \item{A}{
    total association (linear and nonlinear) between Y and the group X. 
  }
  \item{rA}{
    the residual association (the association left in the residuals
    after the linear part has been regressed out of Y). 
  }
  \item{nl1}{
    A - Rsq, the nonlinear part of the association. 
  }
  \item{nl2}{
    (A - Rsq) / A, the nonlinear proportion of the association. 
  }
  \item{nl3}{
    (A - Rsq) / (1 - Rsq), the proportion of total variance that
    is not explained by a linear model but is explained by A.
  }
}
\references{
coming soon
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
\code{\link{agram}}
}

\examples{
    X1 = runif(1000)
    X2 = runif(1000)
    Y = sin(0.5*pi*X1) + sin(0.5*pi*X2) + rnorm(1000)*0.000001
    ma.nl(Y,cbind(X1,X2))
    #
    # in the case of bivariate associations all these measures
    # are symmetric apart from rA, the residual association
    X = runif(1000)
    Y = sin(0.5*pi*X) + rnorm(1000)*0.01
    ma.nl(Y,X)$rA
    ma.nl(X,Y)$rA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
