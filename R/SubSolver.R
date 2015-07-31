#' reformulate sampling data
#'
#' @param matZ randomly generated systematic risk factors, matrix in (n,d) dimensional space
#' @param veceps randomly generated idiosyncratic risk factors, vector in m dimension
#' @param rownum index of obiligor, postive integer
#' @param n number of scenarios, postive integer
#' @return reformulated matrix (Z, eps_k) and x_k.e
Local_data <- function(matZ, veceps, rownum, n){
  # bufA <- cbind(matZ, rep(veceps[rownum], n))
  bufA <- cbind(matZ, rnorm(n,0,1))
  bufb <- rep(veceps[rownum], n)
  reslist <- list(bufA, bufb)
  return(reslist)
}

#' SDP solver for $norm(Ax-b)$
#'
#' @param matA matrix A, matrix in (n,d+1) dimensional space
#' @param recb vector b, vector in n dimension
#' @return optimal matrix X and it's rank
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

#' recovery vector from optimal matrix solution by SDP solver
#'
#' @param matX matrix X, matrix in (d+2,d+2) dimensional space
#' @return approximately recoveried x
sdprecovery_x <- function(matX){
  x <- sign(matX[1,2:(d+2)])*(diag(matX)^(1/2))[2:(d+2)]
  return(x)
}


#' Subfunction for Projection method
#'
#' @param matA matrix A, matrix in (n,d+1) dimensional space
#' @param recb vector b, vector in n dimension
#' @param lambda lagrange multiplier
#' @return list of x(lambda) and y(lambda)
xy_lambda <- function(matA, vecb, lambda){
  dimA <- dim(matA)
  bufmatA <- t(matA) %*% matA + lambda * diag(dimA[2])
  bufvecb <- t(matA) %*% vecb
  x_lambda <- solve(bufmatA, bufvecb)
  y_lambda <- solve(bufmatA, x_lambda)
  xy_lambda <- list(x_lambda, y_lambda)
  return(xy_lambda)
}

#' projection method for soloving $min norm(Ax-b)$
#'
#' @param matA matrix A, matrix in (n,d+1) dimensional space
#' @param recb vector b, vector in n dimension
#' @param lambda_init = 0, set initial lagrange multiplier
#' @param TOL = 1e-4, set termination tolerance
#' @param MATITER = 1e+3, set iteration limitation
#' @return optimal minimum
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

#' find optimal risk loading coefficients
#'
#' @param Zd sampling loading risk vectors, matrix in (n,d) dimensional space
#' @param epsm sampling loading, vector in n dimension
#' @param method = 0, using projection method, method = 1, using SDP method
#' @return list of loading matrix A and loading vector b
#' @export
Loading_coeff <- function(Zd, epsm, n, method = 0){
  loading_matA <- matrix(rep(0,m*d), nrow = m, ncol = d)
  loading_vecb <- rep(0,m)
  for(k in 1:m){
    reslist <- Local_data(Zd, epsm, k, n)
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
