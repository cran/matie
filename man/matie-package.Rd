\name{matie-package}
\alias{matie-package}
\alias{matie}
\docType{package}
\title{
Measuring Association and Testing Independence Efficiently
}
\description{
A non-parametric measure of association between variables. 
The association score, \eqn{A}, ranges from 0 (when the variables are independent) 
to 1 (when they are perfectly associated). 
\eqn{A} is a kind of \eqn{R^2} estimate, and can be thought of as 
the proportion of variance in one variable explained by another 
(or explained by a number of other variables - 
\eqn{A} works for multivariate associations as well).

\code{matie} computes \eqn{A} by estimating a generalized \eqn{R^2}, 
which is computed from the ratio of the likelihood 
of an alternative model (allowing dependence between variables) 
over the likelihood of a null model (that forces the variables to be independent).
See \url{http://arxiv.org/abs/1303.1828} for details.

The variables should be continuous. 
\eqn{A} is calculated from the ranks of the values,
and ties are currently broken randomly. Ordinal data can be handled, but the
properties of \eqn{A} have not been thoroughly investigated for non-continuous 
data with many ties. 
\eqn{A} is not currently implemented for category-valued data.
}
\details{
\tabular{ll}{
Package: \tab matie\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-05-17\cr
License: \tab GPL-3\cr
}
}
\author{
Ben Murrell, Dan Murrell, Hugh Murrell. 
Maintainer: Hugh Murrell <hugh.murrell@gmail.com>
}
\references{
Discovering general multidimensional associations, \url{http://arxiv.org/abs/1303.1828}
}

\keyword{ association }
\keyword{ likelihood }
\keyword{ distribution }
\keyword{ density }
\concept{ maximum likelihood }

\examples{
    d <- sbd(NULL,n=400,Rsq=0.9)
    aScore <- ma(d)
    pValue <- ma.test(d,aScore)
}
