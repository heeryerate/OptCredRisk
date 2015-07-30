########################################################################################
#                                                                                      #
#                                      Packages                                        #
#                                                                                      #
########################################################################################
if(TRUE){
  library(ggplot2)
  require(MASS)
  require(Matrix)
}

########################################################################################
#                                                                                      #
#                                Numerical Experiment                                  #
#                                                                                      #
########################################################################################
if(TRUE){
  NumericalExamples <- function (pk, ck, A, b, x, y, DriftZ = "yes"){
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
}

########################################################################################
#                                                                                      #
#                              Golbal Data & Test Data                                 #
#                                                                                      #
########################################################################################
if(TRUE){
  input <- c(15, 10, 2, 500, 5000)
  assign("m", input[1], envir = .GlobalEnv)        # Number of obligors
  assign("d", input[2], envir = .GlobalEnv)         # Number of factors
  assign("n", input[3]*(d+1), envir = .GlobalEnv)    # Number of samplings
  assign("IS_Sampling_Num", input[4], envir = .GlobalEnv)    # Number of samplings
  assign("PlainMC_Sampling_Num", input[5], envir = .GlobalEnv)    # Number of samplings

  ck <- ceiling((5/m)*(1:m))^2
  # ck <- rep(1,m)
  Zd <- mvrnorm(n, rep(0,d), diag(d))         # Sampling systematic risk factor

#   x <- 100.0
#   y <- seq(20,800,10)
#   ind <- length(y)-5
  x <- 16
  y <- seq(8,140,5)
  ind <- length(y)-5

  lossdata <- list()
  for(i in 1:8){
    if(i <= 4){
        Zd <- mvrnorm(n, rep(0,d), diag(d))         # Sampling systematic risk factor
        epsm <- rnorm(m, mean = 0, sd = 1)
        reslist <- Loading_coeff(Zd, epsm, method = 0)
        A <- reslist[[1]]
        b <- reslist[[2]]

        corrplot(A %*% t(A), method = "circle") #plot matrix

        print("Regression Coefficient Generated")

        lossdata[[i]] <- Numerical_example(pk, ck, A, b, x, y)
        print("Tail probabilities Generated")
    }else{
      A <- matrix(runif(m*d, min=-1, max=1), nrow=m, ncol=d, byrow=T)
      b <- sign(runif(m,-1,1)) * sqrt(1-rowSums(A^2))
      corrplot(A %*% t(A), method = "circle") #plot matrix

      print("Random Coefficient Generated")

      lossdata[[i]] <- Numerical_example(pk, ck, A, b, x, y)
      print("Tail probabilities Generated")
    }
  }
}
