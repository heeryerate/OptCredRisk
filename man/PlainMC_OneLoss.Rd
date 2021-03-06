% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/MC.R
\name{PlainMC_OneLoss}
\alias{PlainMC_OneLoss}
\title{single loss probability $P(L>x)$ by MC simulation}
\usage{
PlainMC_OneLoss(PlainMC_Sampling_Num, m, d, pk, ck, A, b, x)
}
\arguments{
\item{PlainMC_Sampling_Num}{scenarios for MC simiulation, positive integer}

\item{m}{number of obligors, positive integer}

\item{d}{number of risk factors, postive integer}

\item{pk}{default probability of obligors, vector in m dimension}

\item{ck}{expoure of obligors, vector in m dimenstion}

\item{A}{risk loading matrix, matrix in (m, d) dimension}

\item{b}{risk loading vector, vector in m dimension}

\item{x}{threshold of portfolio, postive real value}
}
\value{
average $L$, $P(L>X)$ and 95 percents confidence parameter $sigma$
}
\description{
single loss probability $P(L>x)$ by MC simulation
}

