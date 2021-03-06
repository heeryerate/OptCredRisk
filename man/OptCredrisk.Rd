% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/NumExas.R
\name{OptCredrisk}
\alias{OptCredrisk}
\title{Numerical experiments in OptCredRisk}
\usage{
OptCredrisk(PlainMC_Sampling_Num, IS_Sampling_Num, n, d, m, pk, ck, x, y)
}
\arguments{
\item{PlainMC_Sampling_Num}{scenarios for MC simiulation, positive integer}

\item{IS_Sampling_Num}{scenarios for IS simiulation, positive integer}

\item{n}{number of scenarios, postive integer}

\item{d}{number of risk factors, postive integer}

\item{m}{number of obligors, positive integer}

\item{pk}{default probability of obligors, vector in m dimension}

\item{ck}{expoure of obligors, vector in m dimenstion}

\item{x}{threshold of portfolio, postive real value}

\item{y}{set of threshold of portfolio, postive vector}

\item{A}{risk loading matrix, matrix in (m, d) dimension}

\item{b}{risk loading vector, vector in m dimension}

\item{DriftZ}{flag indicator on whether we need shift mean of risk factors}
}
\value{
Large Loss probability for both IS and MC in optimal coefficients and random coefficients
}
\description{
Numerical experiments in OptCredRisk
}

