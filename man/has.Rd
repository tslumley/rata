\name{has}
\alias{has}
\alias{hasonly}
\alias{hasall}
\alias{hasany}
\alias{wide}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify predictor terms from multiple-response objects
}
\description{
Create a logical vector as a regression predictor that specifies whether a multiple response has particular levels
}
\usage{
has(mr_object, value)
hasonly(mr_object, value)
hasany(mr_object, values)
hasall(mr_object, values)
wide(mr_object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mr_object}{
A multiple-response object
}
  \item{value,values}{
A character vector with one or more levels of the object
}
}

\value{
A logical vector for the 'has' functions; a logical matrix for \code{wide}
}

\seealso{
\code{\link{long}},\code{\link{mrglm}}
The \code{rimu} package.
}
\examples{
data(nzbirds, package="rimu")
birds<-as.mr(nzbirds)
birds
has(birds, "kea")
hasany(birds, c("kea","tui"))
hasall(birds, c("kea", "tui"))

}
