\name{vetinfo}
\alias{vetinfo}
\docType{data}
\title{
Data on sources of veterinary information, from Loughin and Scherer
}
\description{
 Livestock farmers in Kansas were classified as to to their educational level (high school or less, vocational school, 2-year college, 4- other), and they were asked, "What are your primary sources of veterinary info? Options were (A) professional consultant, (B) veterinarian, (C) state or service, (D) magazines, and (E) feed companies and reps.}
\usage{data("vetinfo")}
\format{
  A data frame with 262 observations on the following 2 variables.
  \describe{
    \item{\code{education}}{a numeric}
    \item{\code{sources}}{a matrix with 5 columns}
  }
}

\source{
Loughin TM, Scherer PN (1998) "Testing for Association in Contingency Tables with Multiple Column Responses" Biometrics 54: 630-637
}
\references{
Also analysed by

Agresti A and Liu I-M (1999) "Modeling a Categorical Variable Allowing Arbitrarily Many Category Choices" Biometrics 55:936-943

Decady YJ, Thomas DR (2000) "A Simple Test of Association for Contingency Tables with Multiple Column Responses" Biometrics 56: 893-896
}
\examples{
data(vetinfo)
with(vetinfo, mtable(education, sources))

logistic<-mrglm(present(sources)~value(sources)+education, family=binomial, data=vetinfo)
summary(logistic)
anova(logistic)


loglin<-mrloglin(~sources+education, family=binomial, data=vetinfo)
summary(loglin)
anova(loglin, update(loglin, .~.*education))

}
\keyword{datasets}
