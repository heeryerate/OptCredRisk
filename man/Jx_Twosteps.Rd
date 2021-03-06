% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/SuppFuns.R
\name{Jx_Twosteps}
\alias{Jx_Twosteps}
\title{Approximation upper bound function J(x) for Two-steps Importance Sampling}
\usage{
Jx_Twosteps(Z, d, m, ck, pk, A, b, x)
}
\arguments{
\item{Z}{systematic risk factor, vector in d dimension}

\item{d}{number of risk factors, postive integer}

\item{m}{number of obligors, positive integer}

\item{ck}{expoure of obligors, vector in m dimenstion}

\item{pk}{default probability of obligors, vector in m dimension}

\item{A}{risk loading matrix, matrix in (m, d) dimension}

\item{b}{risk loading vector, vector in m dimension}

\item{x}{threshold of portfolio, postive real value}
}
\value{
function value of J(x)
}
\description{
Approximation upper bound function J(x) for Two-steps Importance Sampling
}

