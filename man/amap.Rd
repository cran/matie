\name{amap}
\alias{amap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap and dendrogram from associations between all pairs of variables. 
}
\description{
Clusters variables in a datafrome using dist = 1-A 
and plots a heatmap with cluster dendrogram attached.
}
\usage{
amap(dataSet, palette, corAdjusted, method, title, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
  the n x m data frame representing n samples of the m variables.
  }
  \item{palette}{
  colour palette to use for the heatmap.
  }
  \item{corAdjusted}{
  a flag which if set causes bi-variate A scores to be mutiplied
  by the sign of the same bi-variate correlation. This gives a
  direction to the association. defaults to \code{FALSE}.
  }
  \item{method}{
  the method to use for correlation association as performed by \code{cor}, 
  one of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}. The default
  is \code{"spearman"}.
  }
  \item{title}{
  a main title for the heatmap plot.
  }
  \item{...}{
  passed on to graphics routines ...
  }
}
\details{
  uses \code{heatmap.2} from the \code{gplots} package
  to display a heatmap and associated dendorgram with
  bi-variate scores provided by matie's A value.
  If required these scores can be expanded into the
  range \code{[-1,1]} through multiplication by the
  sign of the corresponding bi-variate correlation score.
}
\value{
  Returns nothing.
}
\references{
This code was adapted from Moignard et al.
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
\code{\link{tap}}
}

\examples{
    data(baseballData)
    amap(baseballData,corAdjusted=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ baseballData }
\keyword{ association }
