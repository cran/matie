\name{matie-package}
\alias{matie-package}
\alias{matie}
\docType{package}
\title{
Measuring Association and Testing Independence Efficiently
}
\description{
Uses a maximum likelihood density ratio to estimate association between pairs of variables from a data set.
}
\details{
\tabular{ll}{
Package: \tab matie\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-17\cr
License: \tab GPL-3\cr
}
}
\author{
Ben Murrell, Dan Murrell, Hugh Murrell. 
Maintainer: Hugh Murrell <murrellhc@gmail.com>
}
\references{
coming soon
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
