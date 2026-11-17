\name{as.record}
\alias{as.record}
\title{Converts individual capture histories or presence/absence ascertainment records into the categorical input format required by CARE1.}
\usage{as.record(x)}
\description{ Converts individual capture histories or presence/absence ascertainment records into the categorical input format required by CARE1.}
\value{A vector containing the number of individuals in various categories.}
\arguments{\item{x}{A matrix containing the observed individual capture histories or presence/absence ascertainment records.}}
\details{
Each row of the input matrix represents an individual's capture history or presence/absence ascertainment record across multiple sampling occasions or lists. Each matrix element must be 0 (absence or non-capture) or 1 (presence or capture). \cr

For example, consider the following matrix containing the presence/absence records of 10 individuals across three lists. The function classifies these 10 individuals into seven capture-history categories: (001), (010), (011), (100), (101), (110), and (111), and computes the number of individuals in each category. In this example, three individuals were recorded in lists 1 and 3 but not in list 2 (category 101), while no individuals were recorded in all three lists (category 111). Similar counts can be obtained for other categories. \cr

Data for 10 individuals: \cr

\tabular{rlll}{
[1,]  \tab 0 \tab 1 \tab 0 \cr
[2,]  \tab 1 \tab 0 \tab 1 \cr
[3,]  \tab 0 \tab 1 \tab 0 \cr
[4,]  \tab 0 \tab 1 \tab 1 \cr
[5,]  \tab 1 \tab 0 \tab 1 \cr
[6,]  \tab 1 \tab 0 \tab 1 \cr
[7,]  \tab 0 \tab 0 \tab 1 \cr
[8,]  \tab 0 \tab 1 \tab 0 \cr
[9,]  \tab 0 \tab 1 \tab 0 \cr
[10,] \tab 1 \tab 1 \tab 0 \cr
}

The converted data are: \cr

\tabular{rrrrrrr}{
001 \tab 010 \tab 011 \tab 100 \tab 101 \tab 110 \tab 111 \cr
1   \tab 4   \tab 1   \tab 0   \tab 3   \tab 1   \tab 0   \cr
}

For four-list data, individual records are classified into 15 capture-history categories; for five-list data, they are classified into 31 categories.
}
\examples{
# Classify the capture histories of 100 individuals across three lists 
# into seven category data
set.seed(123)
x <- matrix(sample(0:1, 300, replace = TRUE), ncol = 3)
as.record(x)}
