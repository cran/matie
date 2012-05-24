\name{shpd}
\alias{shpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates Sample HyperPlane Data.
}
\description{
  This function generates a sample hyperplane data set. 
}
\usage{
shpd(n,m=2,Rsq=0.7,Ri=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    number of sample points to generate
  }
  \item{m}{
    a hyperplane of dimension m-1 will be generated in m space
  }
  \item{Rsq}{
    the desired coefficient of determination for the hyperplane,
    indicates how far data points will vary from the plane.
  }
  \item{Ri}{
    correlation coefficient for "independent"" variables within
    the hyperplane
  }
}
\details{
  Scatters data around the hyperplane: Vm = sum(V1...Vm-1) + noise where
  noise is adjusted so that the final dataset has coefficient of determination equal to Rsq.
  Variables V1..Vm-1 are assumed to be independent but the user can make them dependent
  on each other by setting the Vi parameter away from zero.
}
\value{
  Returns a n x m data set
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
}

\examples{
    d <- shpd(500,3,Rsq=0.8,Ri=0.01)
    ma(d)$A
    ma(d,partition=list(1,2))$A
    ma(d,partition=list(3,1))$A
    ma(d,partition=list(3,2))$A
    ma(d,partition=list(3,c(1,2)))$A
    cor(d)^2
    # if m=3 and if you have rgl 
    # you can view the data set in 3D 
    # library("rgl")  
    # plot3d(d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{coefficient of determination}
