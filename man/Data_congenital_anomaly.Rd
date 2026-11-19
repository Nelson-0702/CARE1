\name{Data_congenital_anomaly}
\alias{Data_congenital_anomaly}
\docType{data}
\title{A five-list dataset}
\description{
A five-list dataset on congenital anomalies in infants.
}
\usage{
data(Data_congenital_anomaly)
}
\format{
A vector of length 31 containing the numbers of cases in the 31 presence/absence ascertainment categories required by CARE1.
}
\details{
The data were originally reported by Fienberg (1972) and subsequently analyzed by Wittes et al. (1974). They comprise five incomplete lists of infants born with a specific congenital anomaly in Massachusetts, USA, during a specified time period. \cr

The five sources of ascertainment were:
\enumerate{
\item Obstetric records (183 cases).
\item Other hospital records (215 cases).
\item Records maintained by the State Department of Public Health (36 cases).
\item Records maintained by the State Department of Mental Health (263 cases).
\item Records maintained by special schools (252 cases).
}

A total of 537 distinct cases were identified across the five lists.
}
\references{
Fienberg, S. E. (1972). The multiple recapture census for closed populations and incomplete 2k contingency tables. \emph{Biometrika}, 59, 591-603. \cr

Wittes, J. T., Colton, T., and Sidel, V. W. (1974). Capture-recapture methods for assessing the completeness of case ascertainment when using multiple information sources. \emph{Journal of Chronic Diseases}, 27, 25--36.
}
\examples{
data(Data_congenital_anomaly)
estN(Data_congenital_anomaly)
}
\keyword{datasets}
