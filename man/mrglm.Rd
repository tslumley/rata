\name{mrglm}
\alias{mrglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit marginal generalised linear models to multiple-response data
}
\description{
This function fits generalised linear models using \code{\link{glm}}, after expanding the data to represent multiple-response categories in either the outcome or predictors. 
}
\usage{
mrglm(formula, family = gaussian, data, weights, subset, na.action, start = NULL, etastart, mustart, offset, control = list(...), model = TRUE, x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
As in \code{\link{glm}} except for some special functions representing multiple-response objects, see Details. 
}
  \item{family}{
As in \code{\link{glm}}}
  \item{data}{
A data frame containing all the variables in the formula. 
}
  \item{weights}{
As in \code{\link{glm}}}
  \item{subset}{
As in \code{\link{glm}}
}
  \item{na.action}{
As in \code{\link{glm}}}
  \item{start}{
As in \code{\link{glm}}}
  \item{etastart}{
As in \code{\link{glm}}}
  \item{mustart}{
As in \code{\link{glm}}}
  \item{offset}{
As in \code{\link{glm}}}
  \item{control}{
As in \code{\link{glm}}}
  \item{model}{
As in \code{\link{glm}}}
  \item{x}{
As in \code{\link{glm}}}
  \item{y}{
As in \code{\link{glm}}}
  \item{singular.ok}{
As in \code{\link{glm}}}
  \item{contrasts}{
As in \code{\link{glm}}}
  \item{\dots}{
As in \code{\link{glm}}}
}
\details{
In addition to the functions for creating binary variables (see \code{\link{has}}), there are three special functions for expanding the data.  The two functions \code{present} and \code{values} are used when the outcome variable is a multiple-response variable. When these are used, \code{mrglm} creates a record for every possible level of the variable, with \code{present()} giving binary membership indicators and \code{values()} giving the corresponding level of the factor.   The function \code{each} is used when a predictor variable has multiple responses. When \code{each} is used, \code{mrglm} expands the data to give a row for each level that is actually present and \code{each()} returns a factor with those levels. 

}

\note{
The special functions \code{present}, \code{values}, and \code{each} should never be preceded by a namespace label (eg, do not use \code{rata::each()}) because the \code{terms} function recognises them by name.
}

\value{
Object of class \code{mrglm} inheriting from \code{glm}
}
\references{
Agresti A and Liu I-M (1999) "Modeling a Categorical Variable Allowing Arbitrarily Many Category Choices" Biometrics 55:936-943
}


\seealso{
The vignette "Working models for multiple response data"
}
\examples{
data(vetinfo)
mrglm(present(sources)~value(sources),data=vetinfo, family=binomial)
mrglm(present(sources)~value(sources)+education,data=vetinfo, family=binomial)

}

