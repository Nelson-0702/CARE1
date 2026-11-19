
\name{estN}
\alias{estN}
\title{Population Size Estimation Using Capture-Recapture Models}

\usage{
estN(z, se = TRUE, nboot = 100, pair = TRUE)
}

\description{
Estimates population size using capture-recapture models and prints
three sets of results, as specified in the Value section. Population
size estimates based on each pair of lists and the sample coverage
approach are provided, together with standard errors, 95\% confidence
intervals, and parameter estimates, as applicable.
}

\arguments{
  \item{z}{
    A vector in the categorical input format required by CARE1,
    or an output object returned by \code{as.record}.
  }

  \item{se}{
    A logical variable indicating whether to compute the bootstrap
    standard errors of the population size estimates based on
    the bootstrap method and the corresponding 95\% confidence
    intervals. Default is \code{TRUE}.

    This argument applies only to the sample coverage approach
    and does not affect the standard errors and confidence intervals
    of the Chapman estimates.
  }

  \item{nboot}{
    Number of bootstrap replications, an integer greater than one.
    Default is 100.
  }

  \item{pair}{
    A logical variable indicating whether to compute and print
    population size estimates based on each pair of lists.
    Default is \code{TRUE}.
  }
}

\details{
Three population size estimators are provided under the
sample coverage framework:

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
Prints three sets of results:

\enumerate{
  \item The number of identified cases in each list.

  \item Population size estimates based on each pair of lists,
  including the Petersen and Chapman estimates, the standard
  errors of the Chapman estimates, and the corresponding
  95\% confidence intervals. This set is omitted when
  \code{pair = FALSE}.

  \item Population size estimates based on the sample coverage
  approach, including bootstrap standard errors, 95\%
  confidence intervals, and parameter estimates.
}

The following statistics are reported for the sample coverage approach:

\describe{
  \item{M}{
    Number of distinct individuals ascertained in at least one list.
  }

  \item{D}{
    Average number of individuals in the union of all lists
    except one.
  }

  \item{SC}{
    Estimated sample coverage.
  }

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

  \item{u1, u2, ...}{
    Estimated mean ascertainment probabilities depending
    on the population size estimate.
  }

  \item{r12, r13, ...}{
    Estimated coefficients of covariation (CCV) depending
    on the population size estimate. For three-list data,
    the third-order coefficient \code{r123} is also provided.
  }
}

When \code{se = FALSE}, the \code{se}, \code{cil}, and
\code{ciu} columns for the sample coverage approach contain
\code{NA}. The standard errors and confidence intervals of
the Chapman estimates are still computed when
\code{pair = TRUE}.

The function also prints definitions of the reported statistics
and a rough guideline for estimator selection.

In addition to printing the results, the function invisibly
returns a list containing four components:

\describe{
  \item{\code{list_counts}}{
    A one-row matrix containing the number of identified
    cases in each list.
  }

  \item{\code{pair_estimates}}{
    A matrix containing the Petersen and Chapman estimates,
    the standard errors of the Chapman estimates, and the
    corresponding 95\% confidence intervals for each pair
    of lists. Returns \code{NULL} when \code{pair = FALSE}.
  }

  \item{\code{population_estimates}}{
    A matrix with rows \code{Nhat-0}, \code{Nhat}, and
    \code{Nhat-1}, and columns \code{M}, \code{D},
    \code{SC}, \code{est}, \code{se}, \code{cil},
    and \code{ciu}.
  }

  \item{\code{parameters}}{
    A matrix containing the estimated mean ascertainment
    probabilities and coefficients of covariation
    corresponding to each population size estimate.
  }
}
}

\note{
For population size estimates based on each pair of lists:

\enumerate{
  \item Refer to Seber (1982, pages 59 and 60) for the
  Petersen and Chapman estimators and the standard error
  formula for the Chapman estimator.

  \item A log-transformation is used to construct the
  95\% confidence intervals for the Chapman estimates.
  Refer to Chao (1987) for the construction of the
  confidence interval.
}

For the sample coverage approach:

As a rough guideline, sample coverage is considered sufficiently
high when the estimated sample coverage (\code{SC} in the output)
is at least 55\% and the estimated bootstrap standard error of
\code{Nhat} is less than one-third of the population size estimate.
Otherwise, \code{Nhat-1} is recommended.

When \code{se = FALSE}, the bootstrap standard errors are not
computed, and the standard error criterion cannot be assessed.
}

\references{
Chao, A. (1987). Estimating the population size for
capture-recapture data with unequal catchability.
\emph{Biometrics}, 43: 783--791.

Chao, A. (2015). Capture-recapture for human populations.
\emph{Wiley StatsRef: Statistics Reference Online}.
W. Piegorsch, Ed. (John Wiley & Sons), 1--16.
\doi{10.1002/9781118445112.stat04855.pub2}.

Chao, A., Tsay, P. K., Lin, S. H., Shau, W. Y., and Chao, D. Y.
(2001). The applications of capture-recapture models to
epidemiological data. \emph{Statistics in Medicine},
20, 3123--3157. \doi{10.1002/sim.996}.

Seber, G. A. F. (1982).
\emph{The Estimation of Animal Abundance} (2nd ed.).
Griffin, London.
}

\examples{
# Three-list data
data(Data_HAV)
estN(Data_HAV, se = FALSE)

# Four-list data
data(Data_diabetes)
estN(Data_diabetes, se = FALSE)

# Five-list data
data(Data_congenital_anomaly)
estN(Data_congenital_anomaly, se = FALSE)

# Without pairwise population size estimates
estN(Data_HAV, se = FALSE, pair = FALSE)

\donttest{
# With bootstrap standard errors and confidence intervals
set.seed(123)
estN(Data_HAV, se = TRUE, nboot = 100, pair = TRUE)

# With bootstrap but without pairwise estimates
estN(Data_HAV, se = TRUE, nboot = 100, pair = FALSE)
}
}
