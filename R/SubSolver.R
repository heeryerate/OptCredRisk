########################################################################################
#                                                                                      #        
#                                      Packages                                        #                  
#                                                                                      #
########################################################################################
if(TRUE){
  require(MASS)
  require(Matrix)
}

########################################################################################
#                                                                                      #        
#                                   Data Generated (1)                                 #                  
#                                                                                      #
########################################################################################

if(TRUE){
  # epsm <- rep(0,m)                            # Initial idiosyncratic risk factor
  pk <- rep(0,m)                           # Initial default probability for each obligor
  xk <- rep(0,m)                              # Initial default threshold for each obligor
  
  for (i in 1:m){
    # epsm[i] <- rnorm(1, mean = 0, sd = 1)
    pk[i] <- 0.01*(1+sin((16*pi/m)*i))
    # pk[i] <- 0.01
    xk[i] <- qnorm(1-pk[i])
  }
}

########################################################################################
#                                                                                      #        
#                                  Data Generated (2)                                  #                  
#                                                                                      #
########################################################################################

if(FALSE){
  
  epsm <- rep(0,m)                            # Initial idiosyncratic risk factor
  pk <- rep(0,m)                           # Initial default probability for each obligor
  xk <- rep(0,m)                              # Initial default threshold for each obligor
  
  Timeweibull <- function(alpha, theta, lambda){
    res <- ((-1/alpha)*(log(1-0.99^(1/lambda))))^(1/theta)
    return (res)
  }
  
  defaultprob <- function(t,alpha, theta, lambda){
    pk <- (1-exp(-alpha*(t^theta)))^lambda
    return (pk)
  }
  
  epsk <- function(t,alpha, theta, lambda){
    epsm <- t/Timeweibull(alpha, theta, lambda)
    return (epsm)
  }
  
  alpham <- runif(m,0.01,0.5)
  thetam <- runif(m, 0.1,15)
  lambdam <- runif(m,1,10)
  
  for (i in 1:m){
    epsm[i] <- epsk(1,alpham[i], thetam[i], lambdam[i])
    pk[i] <- defaultprob(1.2,alpham[i], thetam[i], lambdam[i])
    xk[i] <- qnorm(1-pk[i])
  }
}

########################################################################################
#                                                                                      #        
#                               Solver Data Generated                                  #                  
#                                                                                      #
########################################################################################

if(TRUE){
  Local_data <- function(matZ, veceps, rownum){
    # bufA <- cbind(matZ, rep(veceps[rownum], n))
    bufA <- cbind(matZ, rnorm(n,0,1))
    bufb <- rep(xk[rownum], n)
    reslist <- list(bufA, bufb)
    return(reslist)
  }
}

########################################################################################
#                                                                                      #        
#                 SDP relaxation Method to Solve Regression Subproblem                 #                  
#                                                                                      #
########################################################################################
if(TRUE){
  sdp_solver <- function(matA, recb){
    obj_mat <- list(-1*rbind(cbind(t(recb)%*%recb,-t(recb)%*% matA),cbind(-t(matA)%*% recb,t(matA)%*% matA)))
    cons_mat <- list(list(as.matrix(Diagonal(x = c(0,rep(1,d+1))))),list(as.matrix(Diagonal(x = c(1,rep(0,d+1))))))
    block_num <- c(1,1)
    block_type <- list(type=c("s"),size=c(d+2))
    res <- Rcsdp::csdp(obj_mat, cons_mat, block_num, block_type)
    opti_X <- matrix(unlist(res$X),d+2,d+2)
    opti_X_rank <- matrixcalc::matrix.rank(opti_X)
    reslist <- list(opti_X, opti_X_rank)
    return (reslist)
  }
  
  sdprecovery_x <- function(matX){
    x <- sign(matX[1,2:(d+2)])*(diag(matX)^(1/2))[2:(d+2)]
    return(x)
  }
}

########################################################################################
#                                                                                      #        
#                   Projection Method to Solve Regression Subproblem                   #                  
#                                                                                      #
########################################################################################
if(TRUE){ 
  # Define function to return x(\lambda) and y(\lambda)
  xy_lambda <- function(matA, vecb, lambda){
    dimA <- dim(matA)
    bufmatA <- t(matA) %*% matA + lambda * diag(dimA[2])  
    bufvecb <- t(matA) %*% vecb
    x_lambda <- solve(bufmatA, bufvecb)
    y_lambda <- solve(bufmatA, x_lambda)
    xy_lambda <- list(x_lambda, y_lambda)
    return(xy_lambda)
  }
  
  Projection_Main <- function(matA, vecb, lambda_init = 0, TOL = 1e-4, MAXITER = 1e+3){
    lambda <- lambda_init
    iter <- 0
    while(TRUE){
      iter <- iter + 1
      xy <- xy_lambda(matA, vecb, lambda)
      x <- unlist(xy[1])
      y <- unlist(xy[2])
      xdotx <- (x %*% x)[1]
      ydoty <- (y %*% y)[1]
      xdoty <- (x %*% y)[1]
      if(abs(norm(x,type = "2") - 1) <= TOL || iter >= MAXITER){
        return (x)
        break
      }
      delta <- xdoty * xdoty + (xdotx - 1) * xdotx * ydoty
      if (delta <= 0){
        lambda <- lambda - xdoty / ydoty
      }else{
        lambda <- lambda + (sqrt(delta) - xdoty) / ydoty
      }
    }
  }
}

########################################################################################
#                                                                                      #        
#                         Loading coefficient generating                               #                  
#                                                                                      #
########################################################################################
if(TRUE){
  Loading_coeff <- function(Zd, epsm, method = 0){
    loading_matA <- matrix(rep(0,m*d), nrow = m, ncol = d)
    loading_vecb <- rep(0,m)
    for(k in 1:m){
      reslist <- Local_data(Zd, epsm, k)
      matA <- reslist[[1]]
      recb <- reslist[[2]]
      
      if(method == 0){
        vecak_bk <- Projection_Main(matA, recb)
        # print(norm((matA %*% vecak_bk - recb), type = "F"))
      }else{
        reslist <- sdp_solver(matA, recb)
        opti_X <- reslist[[1]]
        
        vecak_bk <- sdprecovery_x(opti_X)
        print(norm((matA %*% vecak_bk - recb), type = "F"))
      }
      loading_vecb[k] <- vecak_bk[d+1]
      loading_matA[k,1:d] <- vecak_bk[1:d]
    }
    reslist <- list(loading_matA, loading_vecb)
    return(reslist)
  }
}

########################################################################################
#                                                                                      #        
#                                    TEST EXAMPLE                                      #                  
#                                                                                      #
########################################################################################
if(FALSE){
  Subprob_solver_test <- function(Zd, epsm, rownum = 1){
    reslist <- Local_data(Zd, epsm, rownum)
    matA <- reslist[[1]]
    recb <- reslist[[2]]
    
    reslist <- sdp_solver(matA, recb)
    opti_X <- reslist[[1]]
    
    xsdp <- sdprecovery_x(opti_X)
    print(xsdp)
    
    xproj <- Projection_Main(matA, recb)
    print(xproj)
  }
  
  loading_generate_test <- function(Zd, epsm, method = 0){
    reslist <- Loading_coeff(Zd, epsm, method)
    matA <- reslist[[1]]
    recb <- reslist[[2]]
  
    print(matA)
    print(recb)
  }
}