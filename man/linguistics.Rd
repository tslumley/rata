\name{linguistics}
\alias{linguistics}
\docType{data}
\title{
Evaluating English statements by non-native speakers
}
\description{
Six raters rated 50 non-native English utterances on a 1-9 scale for overall comprehension and also noted whether there was a problem with (a) pronunciation of consonants, (b) pronunciation of vowels, (c) word stress, (d) sentence stress, (e) rhythm, (f) intonation, or (g) rate. 
}
\usage{data("linguistics")}
\format{
  A data frame with 300 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{identifier for observation}
    \item{\code{utterance}}{identifier for utterance}
    \item{\code{rater}}{identifier for rater}
    \item{\code{rating}}{ease of understanding: 1 is low,  9 is high}
    \item{\code{problem}}{\code{mr} object: presence of the 7 problems listed as a-g in the description}
  }
}

\source{
Liu I, Suesse T (2008) "The Analysis of Stratified Multiple Responses" Biometrical Journal 50:135-149
}

\examples{
data(linguistics)
linguistics$rate3 <- with(linguistics, ifelse(rating<3, 1, ifelse(rating<7, 2,3)))
with(linguistics, mtable(problem, rate3))

\donttest{	  
if (getRversion() >= "4.1.0" && require("ggmosaic")) {
   with(linguistics, mr_stack(problem, rate3)) |>	
     ggplot()+geom_mosaic(aes(x=product(rate3,problem), fill=problem))
}
}

mrglm(present(problem)~0+value(problem), data=linguistics, family=binomial())
mrglm(present(problem)~0+value(problem)+factor(rater), data=linguistics, family=binomial())

mrglm(present(problem)~0+value(problem)*rate3+factor(rater), data=linguistics, family=binomial())

mrglm(rating~wide(problem)+factor(rater), data=linguistics)

}
\keyword{datasets}
