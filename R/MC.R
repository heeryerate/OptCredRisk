#' single loss probability $P(L>x)$ by MC simulation
#'
#' @param PlainMC_Sampling_Num scenarios for MC simiulation, positive integer
#' @param m number of obligors, positive integer
#' @param d number of risk factors, postive integer
#' @param pk default probability of obligors, vector in m dimension
#' @param ck expoure of obligors, vector in m dimenstion
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param x threshold of portfolio, postive real value
#' @return average $L$, $P(L>X)$ and 95 percents confidence parameter $sigma$
#' @export
PlainMC_OneLoss <- function (PlainMC_Sampling_Num, m, d, pk, ck, A, b, x){
  Lx <- rep(0,PlainMC_Sampling_Num);
  pL <- rep(0,PlainMC_Sampling_Num)

  for(iter in 1:PlainMC_Sampling_Num){
    Z <- mvrnorm(1, rep(0,d), diag(d))
    epsk <- rnorm (m, mean = 0, sd = 1)

    xK = NULL
    for (k in 1:m){
      xK[k] <- A[k, ] %*% Z + b[k]*epsk[k]
      ############################### t-copula #############################
      #          w <- rnorm(1,mean = 0, sd = 1)
      #          xK[k] <- (A[k, ] %*% Z + b[k]*epsk[k])/w
      ######################################################################
    }

    yk <- (xK > qnorm(1-pk))
    ############################## Alternate MC ############################
    #      buff_uk <- runif(m, min = 0, max = 1)
    #      yk <- (buff_uk <= pk)
    ########################################################################

    Lx[iter] <- ck %*% yk
    pL[iter] <- (Lx[iter] > x)
  }

  L <- sum(Lx)/PlainMC_Sampling_Num
  pLx <- sum(pL)/PlainMC_Sampling_Num
  sigma <- sqrt(sum((pL-pLx)^2)/(PlainMC_Sampling_Num-1))

  out <- c(L, pLx, sigma)

  return (out)
}

#' loss probability distribution $P(L>x)$ by MC simulation
#'
#' @param PlainMC_Sampling_Num scenarios for MC simiulation, positive integer
#' @param m number of obligors, positive integer
#' @param d number of risk factors, postive integer
#' @param pk default probability of obligors, vector in m dimension
#' @param ck expoure of obligors, vector in m dimenstion
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param y set of threshold of portfolio, postive vector
#' @return list of average $L$, $P(L>X)$ and 95 percents confidence interval $p-$ and $p+$
#' @export
PlainMC_LossSet <- function (PlainMC_Sampling_Num, m, d, pk, ck, A, b, y){

  leny <- length(y)
  L <- rep(0,leny)
  pLy <- rep(0,leny)
  pLyminus <- rep(0,leny)
  pLyplus <- rep(0,leny)
  sigma <- rep(0,leny)
  delta <- 0.05

  Zdelta <- qnorm(1-delta/2)

  for(i in 1:leny){
    out <- PlainMC_OneLoss(PlainMC_Sampling_Num, m, d, pk, ck, A, b, y[i])
    L[i] <- out[1]
    pLy[i] <- out[2]
    sigma[i] <- out[3]
  }

  pLyminus <- pLy - Zdelta*sigma/sqrt(PlainMC_Sampling_Num)
  pLyplus <- pLy + Zdelta*sigma/sqrt(PlainMC_Sampling_Num)

  out <- list(L, pLy, pLyminus, pLyplus)

  return (out)
}
