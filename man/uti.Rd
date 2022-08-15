\name{uti}
\alias{uti}
\docType{data}
\title{
Urinary tract infections and contraception choice
}
\description{
These data describe presence of urinary tract infections in 239 female university students and the type(s) of contraception they used.
}
\usage{data("uti")}
\format{
  A data frame with 240 observations on the following 4 variables.
  \describe{
    \item{\code{case}}{identifier}
    \item{\code{age}}{factor with levels "<24" and "24+"}
    \item{\code{uti}}{urinary tract infection: factor with levels "yes" and "no"}
    \item{\code{item}}{\code{mr} object with codes a: oral contraceptive, b: condom, c:lubricated condom, d: spermicide, e:diaphragm}
  }
}

\details{
 Liu and Suesse attribute the data to Bilder and Loughin, who attribute them to a Cytel software manual. The original source appears to be Foxman et al, who describe this as a case-control study carried out at the University of Michigan, but they have 437 observations rather than 239.  The 240th observation was added by Liu and Suesse to avoid having an empty cell. The original case-control study was assessing whether lubricant/spermicide on condoms affected the risk of urinary tract infections. 
}


\source{
Bilder CR, Loughin TM. Testing for conditional multiple marginal independence. Biometrics. 2002 Mar;58(1):200-8
}
\references{
Foxman B, Marsh J, Gillespie B, Rubin N, Koopman JS, Spear S. Condom use and first-time urinary tract infection. Epidemiology. 1997 Nov;8(6):637-41

Liu I, Suesse T (2008) "The Analysis of Stratified Multiple Responses" Biometrical Journal 50:135-149
}
\examples{
data(uti)
with(uti, mtable(interaction(uti,age), item))

\dontrun{
if (getRversion()>="4.1.0" && require("ggmosaic")) {
with(uti, mr_stack(uti,age,item)) |>
   ggplot()+geom_mosaic(aes(x=product(uti,item,age),fill=item))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


with(subset(uti, age=="24+"), mr_stack(uti,item)) |>
    ggplot()+geom_mosaic(aes(x=product(item,uti),fill=item))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
}
}

mrglm(I(uti=="yes")~wide(item)+age, data=uti, family=binomial)

mrglm(I(uti=="yes")~each(item)+age, data=uti, family=binomial)
}
\keyword{datasets}
