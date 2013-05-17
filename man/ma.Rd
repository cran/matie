\name{ma}
\alias{ma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure association
}
\description{
A non-parametric measure of association between variables. 
The association score \eqn{A} ranges from 0 (when the variables are independent) 
to 1 (when they are perfectly associated). 
\eqn{A} is a kind of \eqn{R^2} estimate, 
and can be thought of as the proportion of variance in one variable 
explained by another 
(or explained by a number of other variables - 
\eqn{A} works for multivariate associations as well).
}
\usage{
ma(d,partition,ht,hp,hs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
  the \code{n x m} data frame containing \code{n} observations of \code{m} variables 
  for which the maximal joint/marginal likelihood ratio score is required.
  }
  \item{partition}{
  a list of column indices specifying variable groupings. 
  
  Defaults to \code{list(c(m),c(1:m-1))} where \code{m = ncol(d)} 
  which indicates explaining the
  last variable by means of all the other variables in the data set.
  }
  \item{ht}{
  tangent for the hyperbolic correction, default \code{ht = 43.6978644}. 
  }
  \item{hp}{
  power for the hyperbolic correction, default \code{hp = 0.8120818}.
  }\item{hs}{
  scale for the hyperbolic correction, default \code{hs = 6.0049711}.
  }
}
\details{
  An estimate of association (possibly nonlinear) is computed 
  using a ratio of 
  maximum likelihoods for the \emph{marginal distribution} and 
  maximum weighted likelihoods for the \emph{joint distribution}.
  
  Before the computation is carried out the data is ranked using the
  \code{rwt} function from the \code{matie} package.
  This estimate is usually conservative (ie low) and a small-samples hyperbolic 
  correction is applied by adding an offset, \code{os},
  to the joint likelihood given by:
  
  \if{latex}{
   \eqn{\code{os} = \left( 1 - \frac{1}{1 + \code{A} \times \code{ht}} \right)
   \frac{\code{n}^{\code{hp}}}{ \code{hs} } }
  }
  
  \if{html}{
   \code{os = ( 1 - 1 / (1 + A * ht) ) * ( n ^ (hp) / hs ) }
  }
   
   before the likelihood ratio is re-computed.
}
\value{
  Returns a list of values ...

  \item{A }{a score (including hyperbolic correction) estimating association for the data}
  \item{rawA}{the association score before hyperbolic correction}
  \item{jointKW }{the optimal kernel width for the joint distribution}
  \item{altLL }{the optimal weighted log likelihood for the alternate distribution}
  \item{nullLL }{the optimal log likelihood for the marginal distribution}
  \item{marginalKW }{the optimal kernel width for the marginal distribution}
  \item{weight}{the optimal weight used for the mixture}
  \item{LRstat}{the \code{LR} statistic, required for computing \code{p} values.}
  \item{nRows}{n, the number of complete samples in the data set}
  \item{mCols}{m, the number of variables in the data set}
  \item{partition}{user supplied partition for the variables in the data set}
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
\code{\link{rwt}}
\code{\link{pd}}
\code{\link{sbd}}
\code{\link{shpd}}
\code{\link{std}}
}

\examples{
    # bivariate association
    d <- shpd(n=1000,m=2,Rsq=0.9)
    ma(d)
    #
    # multivariate association (the proportion of variance in "Salary"
    # explained by "Hits" and "Years")
    data(baseballData)
    ma(baseballData,partition=list(11,c(2,7)))$A
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }
\keyword{ correlation }
\keyword{ nonlinear }
\keyword{ distribution }
