########################################################################################
#                                                                                      #        
#                                      Packages                                        #                  
#                                                                                      #
########################################################################################
if(TRUE){
  require(MASS)
}

########################################################################################
#                                                                                      #        
#                           IS Two-step for single loss                                #                  
#                                                                                      #
########################################################################################
if(TRUE){
  IS_TwoSteps_OneLoss_for_x <- function (IS_Sampling_Num, m, d, pk, ck, A, b, x, y, DriftZ){
    if(DriftZ == "yes"){
      Drift_mu <- Max_Jx_TwoSteps(d, m, ck, pk, A, b, x)
    }else{
      Drift_mu <- rep(0,d)
    }
    
    Lx <- rep(0,IS_Sampling_Num)
    pL <- rep(0,IS_Sampling_Num)
    
    for(iter in 1:IS_Sampling_Num){
      Z <- mvrnorm(1, Drift_mu, diag(d))
      
      pkZ = rep(0,m)
      for (i in 1:m){
        ################################ t-copula model #########################
        #          w <- rnorm(1,mean = 0, sd = 1)
        #          pkZ[i] <- pnorm ((A[i, ] %*% Z + qnorm(pk[i])*w)/b[i])
        #########################################################################
        
        if(b[i] > 0 ){
          pkZ[i] <- pnorm ((A[i, ] %*% Z + qnorm(pk[i]))/b[i])
        }else{
          pkZ[i] <- 1 - pnorm ((A[i, ] %*% Z + qnorm(pk[i]))/b[i])
        }
      }
      
      ELZ <- sum(pkZ*ck)
      if(ELZ >= x){
        theta <- 0
      }else{
        theta <- findroot_twisting_psitheta_equation(x, ck, pkZ)
      }
      
      prob <- twisting_probability_psitheta(theta, x, ck, pkZ)
      
      ############################## Alternate MC ############################
      #       xK = NULL
      #       for (k in 1:m){
      #         epsk <- rnorm (1, mean = 0, sd = 1)
      #         xK[k] <- A[k, ] %*% Z + b[k]*epsk
      #       }
      #       yk <- (xK > qnorm(1-prob))
      ########################################################################
      
      buff_uk <- runif(m, min = 0, max = 1)
      yk <- (buff_uk <= prob) 
      
      Lx[iter] <- ck %*% yk
      buff_phiZ <- psitheta(theta, x, ck, pkZ)
      
      pL[iter] <- (Lx[iter] > y)*exp(-theta*Lx[iter] + buff_phiZ)*exp(-Drift_mu %*% Z + Drift_mu %*% Drift_mu /2)
    }
    
    L <- sum(Lx)/IS_Sampling_Num
    pLy <- sum(pL)/IS_Sampling_Num
    sigma <- sqrt(sum((pL-pLy)^2)/(IS_Sampling_Num-1))
    
    out <- c(L, pLy, sigma)
    return (out)
  }
}

########################################################################################
#                                                                                      #        
#                           IS Two-step for multiple loss                              #                  
#                                                                                      #
########################################################################################
if(TRUE){
  IS_TwoSteps_LossSet <- function (IS_Sampling_Num, m, d, pk, ck, A, b, x, y, DriftZ){
    
    leny <- length(y)
    
    L <- rep(0,leny)
    pLy <- rep(0,leny)
    sigmay <- rep(0,leny)
    pLyminus <- rep(0,leny)
    pLyplus <- rep(0,leny)
    delta <- 0.05
    
    Zdelta <- qnorm(1-delta/2)
    
    for (iter in 1:leny){
      out <- IS_TwoSteps_OneLoss_for_x(IS_Sampling_Num, m, d, pk, ck, A, b, x, y[iter], DriftZ)
      L[iter] <- out[1]
      pLy[iter] <- out[2]
      sigmay[iter] <- out[3]
    }
    
    pLyminus <- pLy - Zdelta*sigmay/sqrt(IS_Sampling_Num)
    pLyplus <- pLy + Zdelta*sigmay/sqrt(IS_Sampling_Num)
    
    out <- list(L, pLy, pLyminus, pLyplus)
    
    return (out)
  }
}