\name{estN.stat}
\alias{estN.stat}
\title{Basic Capture-Recapture Statistics}
\usage{
estN.stat(z)
}
\arguments{
\item{z}{A vector in the categorical input format required by CARE1, or an output object returned by \code{as.record}.}
}

\description{
Calculates the number of identified cases in each list together
with the total number of distinct observed cases (M), the average
number of cases in the union of all lists except one (D), and
the estimated sample coverage (SC).
}

\value{
A one-row numeric matrix containing \code{n1}, ..., \code{nK},
\code{M}, \code{D}, and \code{SC}. List counts and M are unrounded;
D and SC are rounded to three decimal places.
}

\examples{
data(Data_HAV)
estN.stat(Data_HAV)
data(Data_diabetes)
estN.stat(Data_diabetes)
data(Data_congenital_anomaly)
estN.stat(Data_congenital_anomaly)
}
