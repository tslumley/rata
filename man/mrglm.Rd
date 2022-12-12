\name{mrglm}
\alias{mrglm}
\alias{anova.mrglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit marginal generalised linear models to multiple-response data
}
\description{
This function fits generalised linear models using \code{\link{glm}}, after expanding the data to represent multiple-response categories in either the outcome or predictors. 
}
\usage{
mrglm(formula, family = gaussian, data, weights, subset, na.action, start = NULL, etastart, mustart, offset, control = list(...),
model = TRUE, x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...)
\method{anova}{mrglm}(object, object2 = NULL, test = c("F", "Chisq"), method = c("LRT", 
    "Wald"), tolerance = 1e-05, ...) 
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
  \item{weights, subset, na.action, start, etastart,mustart,offset, control,model,x,y,singular.ok,contrasts}{
As in \code{\link{glm}}}
\item{object,object2}{Objects of class \code{mrglm} representing nested models to be compared}
\item{method,test}{Rao-Scott working likelihood ratio test or Wald test, and F or chi-squared approximation. F is usually preferable}
\item{tolerance}{Tolerance for checking whether models are actually nested}
\item{\dots}{not used}
}
\details{
In addition to the functions for creating binary variables (see \code{\link{has}}), there are three special functions for expanding the data.  The two functions \code{present} and \code{values} are used when the outcome variable is a multiple-response variable. When these are used, \code{mrglm} creates a record for every possible level of the variable, with \code{present()} giving binary membership indicators and \code{values()} giving the corresponding level of the factor.   The function \code{each} is used when a predictor variable has multiple responses. When \code{each} is used, \code{mrglm} expands the data to give a row for each level that is actually present and \code{each()} returns a factor with those levels.

Only one variable may have \code{present()}/\code{values()} terms. This may be relaxed in future versions if there are good examples of it being useful.

}

\note{
The special functions \code{present}, \code{values}, and \code{each} should never be preceded by a namespace label (eg, do not use \code{rata::each()}) because the \code{terms} function recognises them by name.
}

\value{
Object of class \code{mrglm} inheriting from \code{glm}
}
\references{
Agresti A and Liu I-M (1999) "Modeling a Categorical Variable Allowing Arbitrarily Many Category Choices" Biometrics 55:936-943

Lumley T, Scott A (2014) "Tests for Regression Models Fitted to Survey Data"  Australian and New Zealand Journal of Statistics 56: 1-14
}


\seealso{
The vignette "Working models for multiple response data"

\code{\link{anova.svyglm}}
}
\examples{
data(vetinfo)
(model0<-mrglm(present(sources)~value(sources),data=vetinfo, family=binomial))
(model1<-mrglm(present(sources)~value(sources)+education,data=vetinfo, family=binomial))
model2<-mrglm(present(sources)~value(sources)*education,data=vetinfo, family=binomial)
anova(model1,model0)
anova(model2)


mrglm(as.integer(education)~each(sources),data=vetinfo)
mrglm(as.integer(education)~wide(sources),data=vetinfo)


}

