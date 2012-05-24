\name{fdg}
\alias{fdg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Force Directed Graph from all pairs association.
}
\description{
Constructs and draws a force directed graph using all pairs
association matrix to drive attraction in the graph.
This adjacency matrix can represent linear or non-linear
associations.
}
\usage{
fdg(dataSet, dataName, method, cutoff, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
  a data frame
  }
  \item{dataName}{
  a string naming the dataSet. (set to NULL if unknown)
  }
  \item{method}{
  one of "A" for non-linear associations or "Rsq" for linear associations.
  Defaults to "A".
  }
  \item{cutoff}{
  a real number between 0 and 1, used to remove edges from the graph with
  association < cutoff.
  Defaults to \code{0.1}
  }
  \item{dim}{
  graph layouts can be generated in 2 or 3 dimensions. note that \code{rgl}
  must be installed for viewing 3 dimensional layouts.
  }
}
\details{
  The \code{igraph} package is used to generate, layout and display the force 
  directed graph. The layout employed is \code{\link{layout.kamada.kawai}}.
}
\value{
  Returns nothing.
}
\references{
coming soon
}
\author{
Ben Murrell, Dan Murrell & Hugh Murrell.
}
\note{
All non numeric columns and non complete rows are stripped from the data set
before the adjacency matrix is computed.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ma}}
\code{\link{cor}}
\code{\link{agram}}
}

\examples{
    data(baseballData)
    fdg(baseballData,dataName="baseballData",method="A",cutoff=0.35,dim=2)
    
    # if you have rgl available then try dim=3 
    # fdg(baseballData,dataName="baseballData",method="Rsq",cutoff=0.15,dim=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
