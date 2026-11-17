\name{estN.para}
\alias{estN.para}
\title{Estimates of Ascertainment Probabilities and Covariation Coefficients}
\usage{
estN.para(z, nhat)
}
\description{
Reports estimated mean ascertainment probabilities and coefficients of covariation (CCV)
for a specified population size estimate.
}
\arguments{
  \item{z}{A vector in the categorical input format required by CARE1, or an output object returned by \code{as.record}.}
  \item{nhat}{A single numerical population size estimate, such as one element of \code{estN(z)$estimates}.}
}
\value{
A one-row matrix with columns \code{u1}, \code{u2}, etc. (estimated mean ascertainment
probabilities) and \code{r12}, \code{r13}, etc. (estimated coefficients of covariation).
Three-list data additionally have the column \code{r123}.
}
\references{
Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y. (2001).
The applications of capture-recapture models to epidemiological data.
\emph{Statistics in Medicine}, 20, 3123--3157.
}
\examples{
data(Data_HAV)
fit <- estN(Data_HAV)
estN.para(Data_HAV, fit$estimates["Nhat", "est"])
}
