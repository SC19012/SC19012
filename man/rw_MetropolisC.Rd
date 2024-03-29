% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rw_MetropolisC}
\alias{rw_MetropolisC}
\title{Implement a randomwalk Metropolis sampler for generating the standard Laplace distribution}
\usage{
rw_MetropolisC(sigma, x0, N)
}
\arguments{
\item{sigma}{varience}

\item{x0}{initial value}

\item{N}{times of generating sample}
}
\value{
a random sample of size \code{n}
}
\description{
Implement a randomwalk Metropolis sampler for generating the standard Laplace distribution
}
\examples{
\dontrun{
rw <- rw_MetropolisC(1,25,1000)
}
}
