% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_TraiSIE_sim_cpp.R
\name{get_TraiSIE_sim_cpp}
\alias{get_TraiSIE_sim_cpp}
\title{Simulation function to create simulated data as observed data in ABC.}
\usage{
get_TraiSIE_sim_cpp(parameters, K, replicates)
}
\arguments{
\item{parameters}{A vector for CES rates.}

\item{K}{Carrying capacity, Inf for diverdity-independent models.}

\item{replicates}{The number of replicates(islands) for TraiSIE simulation.}
}
\value{
A list contains simulated islands
}
\description{
Simulation function to create simulated data as observed data in ABC.
}
\author{
Shu Xie
}
