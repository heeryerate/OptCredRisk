#' Approximation upper bound function J(x) for Two-steps Importance Sampling
#'
#' @param Z systematic risk factor, vector in d dimension
#' @param d number of risk factors, postive integer
#' @param m number of obligors, positive integer
#' @param ck expoure of obligors, vector in m dimenstion
#' @param pk default probability of obligors, vector in m dimension
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param x threshold of portfolio, postive real value
#' @return function value of J(x)
Jx_Twosteps <- function(Z, d, m, ck, pk, A, b, x){
  pkZ = NULL
  for (i in 1:m){
    pkZ[i] <- pnorm ((A[i, ] %*% Z + qnorm(pk[i]))/b[i])
  }

  EL_conditional_Z <- sum(pkZ*ck)
  if(EL_conditional_Z >= x){
    theta <- 0
  }else{
    theta <- findroot_twisting_psitheta_equation(x,ck, pkZ)
  }

  out <- (-1)*(psitheta(theta, x, ck, pkZ) - theta*x - Z %*% Z/2)

  return(out)
}

#' Maximum of approximation upper bound function J(x) for Two-steps Importance Sampling
#'
#' @param d number of risk factors, postive integer
#' @param m number of obligors, positive integer
#' @param ck expoure of obligors, vector in m dimenstion
#' @param pk default probability of obligors, vector in m dimension
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param x threshold of portfolio, postive real value
#' @return maximum of J(x)
Max_Jx_TwoSteps <- function(d, m, ck, pk, A, b, x){
  out <- optim(par=rep(0,d), fn=Jx_Twosteps, d, m, ck, pk, A, b, x
               , gr = NULL,
               method = c("CG"),
               lower = -Inf, upper = Inf,
               control = list(), hessian = F)
  return(out$par)
}

#' function value $phi(theta)$ for Two-steps Importance Sampling
#'
#' @param theta exponential twisting parameter, postive real value
#' @param ck expoure of obligors, vector in m dimenstion
#' @param x threshold of portfolio, postive real value
#' @param pkZ default probability conditional on Z, postive real value
#' @return function value $phi(theta)$
psitheta <- function(theta, x, ck, pkZ){
  BigM <- 100
  Innerprod_theta_ck <- theta*ck
  Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM

  buff_pek <- pkZ*exp(Innerprod_theta_ck)

  out <- sum(log(1-pkZ+buff_pek))

  return(out)
}

#' function $phi'(theta) - x$ for Two-steps Importance Sampling
#'
#' @param theta exponential twisting parameter, postive real value
#' @param ck expoure of obligors, vector in m dimenstion
#' @param x threshold of portfolio, postive real value
#' @param pkZ default probability conditional on Z, postive real value
#' @return function value $phi'(theta)-x$
twisting_psitheta_equation <- function(theta, x, ck, pkZ){

  BigM <- 100
  Innerprod_theta_ck <- theta*ck
  Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM

  buff_pek <- pkZ*exp(Innerprod_theta_ck)

  out <- sum(ck*buff_pek/(1-pkZ+buff_pek))-x

  return(out)
}

#' find root of equation $phi'(theta) = x$ for Two-steps Importance Sampling
#'
#' @param x threshold of portfolio, postive real value
#' @param ck expoure of obligors, vector in m dimenstion
#' @param pkZ default probability conditional on Z, postive real value
#' @return root of $phi'(theta) = x$
findroot_twisting_psitheta_equation <- function(x, ck, pkZ){
  theta <- try(uniroot(twisting_psitheta_equation,c(0,100), x, ck, pkZ, lower = 0,
                   upper = 100,
                   extendInt = "yes", check.conv = F,
                   tol = 0.0001, trace = 100, maxiter = 100))
  if (class(theta) == "try-error"){
    cat("Catch one theta error, set theta = 0")
    theta <- data.frame(root = 0)
  }
  return(theta$root)
}

#' twisted probability for Two-steps Importance Sampling
#'
#' @param theta exponential twisting parameter, postive real value
#' @param ck expoure of obligors, vector in m dimenstion
#' @param x threshold of portfolio, postive real value
#' @param pkZ default probability conditional on Z, postive real value
#' @return twisted probability by optimal $theta$
twisting_probability_psitheta <- function(theta, x, ck, pkZ){

  BigM <- 100
  Innerprod_theta_ck <- theta*ck
  Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM

  buff_pek <- pkZ*exp(Innerprod_theta_ck)

  out <- buff_pek/(1-pkZ+buff_pek)

  return(out)
}

#' Truncate a matrix into [-1,1]
#'
#' @param matA an input matrix
#' @return output matrix
trunMat <- function(matA){
  dimA <- dim(matA)
  matA[matA > 1] <- 1
  matA[matA < -1] <- -1
  return(matA)
}

# this is a sublime git test
#  I deleted this line to see what will happen. return(matA)
