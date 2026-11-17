\name{estN.pair}
\alias{estN.pair}
\title{Pairwise estimation}
\usage{estN.pair(z)}
\description{Estimating population size based on any pair of samples}
\value{A matrix with one row per pair of lists and columns \code{Petersen}, \code{Chapman}, \code{se}, \code{cil}, and \code{ciu}.}
\note{
	\enumerate{
		\item Refer to Seber (1982, pages 59 and 60) for the Petersen estimator and the Chapman estimator as well as s.e. formula.
		\item A log-transformation is used to obtain the confidence interval so that the lower limit is always greater than
				the number of ascertained. Refer to Chao (1987, Biometrics,43, 783-791) for the construction of the confidence interval.}
	}

\references{
Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, 43: 783--791. \cr

Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y. (2001). The applications of capture-recapture models to epidemiological data. \emph{Statistics in Medicine}, 20: 3123--3157. \cr

Seber, G. A. F. (1982). \emph{The Estimation of Animal Abundance} (2nd ed.). Griffin, London.
}

\arguments{\item{z}{A vector in the categorical input format required by CARE1, or an output object returned by \code{as.record}.}}
\examples{
data(Data_diabetes)
estN.pair(Data_diabetes)
}
