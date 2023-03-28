\name{mrmultinom}
\alias{mrmultinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multinomial models for multiple-response data
}
\description{
Models a multiple-response outcome variable as a clustered multinomial variable with each multi-response observation as a cluster.  The \code{wide} and \code{each} special terms from \code{\link{mrglm}} can be used if there is also a multiple-response predictor. In contrast to \code{\link{mrloglin}} the predictors do not need to be categorical. 
}
\usage{
mrmultinom(formula, data, parallel = TRUE~1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula with a \code{mr} or \code{vmr} object on the left-hand side
}
  \item{data}{
data frame to look up the variables
}
\item{parallel}{
Specifies which formula terms have different coefficients for different responses. See \code{\link[VGAM]{multinomial}} for more details 
}
\item{\dots}{
Passed to \code{\link[VGAM]{multinomial}} in the VGAM package
}
}

\value{
An object of class \code{c("mrmultinom","svy_vglm")}
}

\note{
Records with no responses for the outcome variable do not contribute to the multinomial model and will be discarded in fitting. 
}

\seealso{
\code{\link{mrloglin}}, \code{\link[VGAM]{multinomial}}
}
\examples{
data(vetinfo)
mrmultinom(sources~as.numeric(education), data=vetinfo)

mrmultinom(sources~education, data=vetinfo)
}
