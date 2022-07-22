\name{mrglm}
\alias{mrglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalised linear models with multiple-response variables
}
\description{
This function fits a generalised linear model to multiple-response data.  The key difference between this and \code{glm} is the \code{\link{long}} function, which is recognised by \code{mrglm} and interpreted to expand the data vertically so that each response gets its own row.  Sandwich variance estimators are then used to get correct standard errors. 
}
\usage{
mrglm(formula, data, family, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Doesn't actually work yet
}
\value{
an object of class \code{mrglm}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


}
