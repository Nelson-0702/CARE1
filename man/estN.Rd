\name{estN}
\alias{estN}
\title{Population Size Estimation Using Capture-Recapture Models}

\usage{
estN(z, se = FALSE, nboot = 200)
}

\description{
Estimates population size using three capture-recapture estimators:
\code{Nhat-0}, \code{Nhat-1}, and \code{Nhat}. All three estimates are
returned together, without requiring a method argument. Optional
bootstrap standard errors and log-transformed 95\% confidence
intervals are provided.
}

\arguments{
  \item{z}{
    A vector in the categorical input format required by CARE1,
    or an output object returned by \code{as.record}.
  }

  \item{se}{
    A logical variable indicating whether to compute the bootstrap
    standard errors of the population size estimates based on
    the bootstrap method. Default is \code{FALSE}.
  }

  \item{nboot}{
    Number of bootstrap replications, an integer greater than one.
    Default is 200.
  }
}

\details{
\describe{
  \item{Nhat-0}{
    Population size estimate under the assumption of independence
    among sources; see Eq. (9) of Chao et al. (2001).
  }

  \item{Nhat-1}{
    One-step population size estimate for cases with low sample
    coverage; see Eq. (17) of Chao et al. (2001), or Eq. (14)
    of Chao (2015) for an explicit formula for three-list data.
  }

  \item{Nhat}{
    Population size estimate for cases with sufficiently high
    sample coverage; see Eq. (16) of Chao et al. (2001), or
    Eq. (13) of Chao (2015) for an explicit formula for
    three-list data.
  }
}
}


\value{
A matrix containing three population size estimates
(\code{Nhat-0}, \code{Nhat-1}, and \code{Nhat})
and the following columns:

\describe{
  \item{est}{
    Population size estimate.
  }

  \item{se}{
    Estimated standard error of the population size estimate
    based on the bootstrap method.
  }

  \item{cil}{
    Lower limit of the 95\% confidence interval,
    calculated using a log-transformation.
  }

  \item{ciu}{
    Upper limit of the 95\% confidence interval,
    calculated using a log-transformation.
  }
}

When \code{se = FALSE}, the \code{se}, \code{cil},
and \code{ciu} columns contain \code{NA}.
}


\note{
As a rough guideline, sample coverage is considered sufficiently
high when the estimated sample coverage (\code{SC} in the output
of \code{estN.stat}) is at least 55\% and the estimated bootstrap
standard error of \code{Nhat} is less than one-third of the
population size estimate. Otherwise, \code{Nhat-1} is recommended.
}


\references{
Chao, A. (2015). Capture-recapture for human populations.
\emph{Wiley StatsRef: Statistics Reference Online}.
W. Piegorsch, Ed. (John Wiley & Sons), 1--16.
\doi{10.1002/9781118445112.stat04855.pub2}.

Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y.
(2001). The applications of capture-recapture models to
epidemiological data. \emph{Statistics in Medicine},
20, 3123--3157. \doi{10.1002/sim.996}.
}

\examples{
data(Data_HAV)
estN(Data_HAV)

data(Data_diabetes)
estN(Data_diabetes)

data(Data_congenital_anomaly)
estN(Data_congenital_anomaly)

\donttest{
set.seed(123)
estN(Data_HAV, se = TRUE, nboot = 200)
}
}
