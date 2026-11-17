\name{Data_HAV}
\alias{Data_HAV}
\docType{data}
\title{A three-list dataset}
\description{
A three-list dataset on hepatitis A virus (HAV) cases.
}
\usage{
data(Data_HAV)
}
\format{
A vector of length 7 containing the numbers of cases in the 7 presence/absence ascertainment categories required by CARE1.
}
\details{
The data were originally analyzed by Chao et al. (2001) and comprise three incomplete lists of hepatitis A virus (HAV) cases at a college in Taiwan. \cr

The three sources of ascertainment are:
\enumerate{
\item P-list: Records based on a serum test taken by the Institute of Preventive Medicine, Department of Health of Taiwan (135 cases).
\item Q-list: Local hospital records reported by the National Quarantine Service (122 cases).
\item E-list: Records collected by epidemiologists (126 cases).
}

After merging the three lists and eliminating duplicate records, a total of 271 distinct cases were identified.
}
\references{
Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y. (2001). The applications of capture-recapture models to epidemiological data. \emph{Statistics in Medicine}, 20: 3123-3157.
}
\examples{
data(Data_HAV)
CARE1.print(Data_HAV)
}
\keyword{datasets}