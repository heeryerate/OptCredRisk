#' Generate Large Loss probability for Numerical experiments in OptCredRisk
#'
#' @param PlainMC_Sampling_Num scenarios for MC simiulation, positive integer
#' @param IS_Sampling_Num scenarios for IS simiulation, positive integer
#' @param m number of obligors, positive integer
#' @param d number of risk factors, postive integer
#' @param pk default probability of obligors, vector in m dimension
#' @param ck expoure of obligors, vector in m dimenstion
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param x threshold of portfolio, postive real value
#' @param y set of threshold of portfolio, postive vector
#' @param DriftZ flag indicator on whether we need shift mean of risk factors
#' @return Large Loss probability for both IS and MC
#' @export
GenerateLossProb <- function (PlainMC_Sampling_Num, IS_Sampling_Num, d, m ,pk, ck, A, b, x, y, DriftZ = "yes"){
  ind <- length(y)-5
  loss <- y[1:ind]

  out <- IS_TwoSteps_LossSet(IS_Sampling_Num, m, d, pk, ck, A, b, x, y, DriftZ)
  Isminus <- log(unlist(out[2]))[1:ind]
  Isplus <- log(unlist(out[4]))[1:ind]
  # plot(y[3:ind],log(unlist(out[2]))[3:ind], type = "l", xlab = "Loss Level", ylab = "Tail probability", col = "red")
  # lines(y[3:ind],log(unlist(out[4]))[3:ind], col = "red")

  out <- PlainMC_LossSet(PlainMC_Sampling_Num, m, d, pk, ck, A, b, y)
  Mcminus <- log(unlist(out[2]))[1:ind]
  Mcplus <- log(unlist(out[4]))[1:ind]
  # lines(y[1:ind],log(unlist(out[2]))[1:ind])
  # lines(y[1:ind],log(unlist(out[4]))[1:ind])

  nameprob <- c(rep("Isminus",ind),rep("Isplus",ind),rep("Mcminus",ind),rep("Mcplus",ind))
  valueprob <- c(Isminus, Isplus, Mcminus, Mcplus)
  colorprob <- c(rep(1,ind),rep(1,ind),rep(2,ind),rep(2,ind))
  typeprob <- c(rep(1,ind),rep(2,ind),rep(3,ind),rep(4,ind))
  lossprob <- data.frame(loss, valueprob, nameprob, colorprob,typeprob)

  return(lossprob)
}


#' Numerical experiments in OptCredRisk
#'
#' @param PlainMC_Sampling_Num scenarios for MC simiulation, positive integer
#' @param IS_Sampling_Num scenarios for IS simiulation, positive integer
#' @param m number of obligors, positive integer
#' @param d number of risk factors, postive integer
#' @param n number of scenarios, postive integer
#' @param pk default probability of obligors, vector in m dimension
#' @param ck expoure of obligors, vector in m dimenstion
#' @param A risk loading matrix, matrix in (m, d) dimension
#' @param b risk loading vector, vector in m dimension
#' @param x threshold of portfolio, postive real value
#' @param y set of threshold of portfolio, postive vector
#' @param DriftZ flag indicator on whether we need shift mean of risk factors
#' @return Large Loss probability for both IS and MC in optimal coefficients and random coefficients
#' @export
OptCredrisk <- function(PlainMC_Sampling_Num, IS_Sampling_Num, n, d, m ,pk, ck, x, y){
  lossdata <- list()
  ind <- length(y)-5
  for(i in 1:8){
    if(i <= 4){
      Zd <- mvrnorm(n, rep(0,d), diag(d))         # Sampling systematic risk factor
      epsm <- rnorm(m, mean = 0, sd = 1)
      reslist <- Loading_coeff(Zd, epsm, n, 0)
      A <- reslist[[1]]
      b <- reslist[[2]]

      corrplot(A %*% t(A), method = "circle") #plot correlation matrix

      print("Regression Coefficient Generated")

      lossdata[[i]] <- GenerateLossProb(PlainMC_Sampling_Num, IS_Sampling_Num, d, m ,pk, ck, A, b, x, y)

      print("Tail probabilities Generated")
    }else{
      A <- matrix(runif(m*d, min=-1, max=1), nrow=m, ncol=d, byrow=T)
      b <- sign(runif(m,-1,1)) * sqrt(1-rowSums(A^2))

      corrplot(A %*% t(A), method = "circle") #plot correlation matrix

      print("Random Coefficient Generated")

      lossdata[[i]] <- GenerateLossProb(PlainMC_Sampling_Num, IS_Sampling_Num, d, m ,pk, ck, A, b, x, y)

      print("Tail probabilities Generated")
    }
  }

  res <- list()
  for(i in 1:8){
    lossprob <- lossdata[[i]]
    lossprob$valueprob[is.infinite(lossprob$valueprob)] <- 0
    lossprob$nameprob <- c(rep("Isminus",ind),rep("Isplus",ind),rep("Mcminus",ind),rep("Mcplus",ind))
    lossprob$Simulation <- factor(lossprob$colorprob, levels = c(1,2), labels = c("IS", "MC"))
    a <- ggplot2::ggplot(lossprob, ggplot2::aes(x = loss, y = valueprob, color = Simulation, shape = nameprob))
    a <- a + ggplot2::scale_x_continuous(name = "Loss Level" , limits = c(10,140)) + ggplot2::scale_y_continuous(name = "log(Tail Probability)", limits = c(-7.5,-2))
    # a <- a + ggplot2::geom_point(size = 3)
    a <- a + ggplot2::geom_line(size = 1.2)
    res[[i]] <- a
    print(a)
  }
  #   gridExtra::grid.arrange(res[[1]],res[[2]],res[[3]],res[[4]])
  #   gridExtra::grid.arrange(res[[5]],res[[6]],res[[7]],res[[8]])
  return(1)
}
