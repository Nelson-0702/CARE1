\name{Data_diabetes}
\alias{Data_diabetes}
\docType{data}
\title{A four-list dataset}
\description{
A four-list dataset on diabetes cases.
}
\usage{
data(Data_diabetes)
}
\format{
A vector of length 15 containing the numbers of cases in the 15 presence/absence ascertainment categories required by CARE1.
}
\details{
The data were originally reported by Bruno et al. (1994) and comprise four incomplete lists of diabetes cases in a community in Italy. \cr

The four sources of ascertainment were:
\enumerate{
\item Diabetic clinic and family physician records (1754 cases).
\item Hospital discharge records (452 cases).
\item Prescription records (1135 cases).
\item Records of purchases of reagent strips and insulin syringes (173 cases).
}

A total of 2069 cases were identified across the four lists.
}

\references{
Bruno, G., LaPorte, R., Merletti, F., Biggeri, A., McCarty, D., and Pagano, G. (1994).
National diabetes programs. Application of capture-recapture to count diabetes?
\emph{Diabetes Care}, 17(6): 548.
}

\examples{
data(Data_diabetes)
estN(Data_diabetes)
}
\keyword{datasets}
