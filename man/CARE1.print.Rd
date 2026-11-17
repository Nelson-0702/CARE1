\name{CARE1.print}
\alias{CARE1.print}
\title{Prints Summary Results for Data in the Categorical Input Format Required by CARE1}
\usage{
CARE1.print(z)
}
\description{
Prints three sets of results, as specified in the Value section.
}
\arguments{
\item{z}{A vector in the categorical input format required by CARE1, or an output object returned by \code{as.record}.}
}
\value{
Prints three sets of results:
\enumerate{
  \item The number of identified cases in each list.
  \item Population size estimates based on each pair of lists.
  \item Population size estimates based on the sample coverage approach, including estimated sample coverage (SC), bootstrap standard errors, 95\% confidence intervals, and parameter estimates.
}

}
\examples{
data(Data_HAV)
\donttest{
CARE1.print(Data_HAV)
}
}
