\name{long}
\alias{long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expand multiple-response data to long form in a model
}
\description{
This function only has any effect in a model formula used in \code{\link{mrglm}}.
It expands the data so there is a row for every response (every combination of responses if there are multiple \code{long} terms in a formula)
}
\usage{
long(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A multiple-response object
} 
}

\value{
Only has an effect in a model formula
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wide}}, \code{\link{mrglm}}
}
\examples{


}

