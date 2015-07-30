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
#                                  J(x) and related                                    #                  
#                                                                                      #
########################################################################################
if(TRUE){
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
  
  Max_Jx_TwoSteps <- function(d, m, ck, pk, A, b, x){
    out <- optim(par=rep(0,d), fn=Jx_Twosteps, d, m, ck, pk, A, b, x
                 , gr = NULL,
                 method = c("CG"),
                 lower = -Inf, upper = Inf,
                 control = list(), hessian = F)
    return(out$par)
  }
}

########################################################################################
#                                                                                      #        
#                                theta(x) and related                                  #                  
#                                                                                      #
########################################################################################
if(TRUE){
  psitheta <- function(theta, x, ck, pkZ){
    BigM <- 100
    Innerprod_theta_ck <- theta*ck 
    Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM
    
    buff_pek <- pkZ*exp(Innerprod_theta_ck)
    
    out <- sum(log(1-pkZ+buff_pek))
    
    return(out)
  }
  
  twisting_psitheta_equation <- function(theta, x, ck, pkZ){
    
    BigM <- 100
    Innerprod_theta_ck <- theta*ck 
    Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM
    
    buff_pek <- pkZ*exp(Innerprod_theta_ck)
    
    out <- sum(ck*buff_pek/(1-pkZ+buff_pek))-x
    
    return(out)
  }
  
  findroot_twisting_psitheta_equation <- function(x, ck, pkZ){
    theta <- uniroot(twisting_psitheta_equation,c(0,100), x, ck, pkZ, lower = 0, 
                     upper = 100, 
                     extendInt = "yes", check.conv = F, 
                     tol = 0.0001, trace = 100, maxiter = 100)
    return(theta$root)
  }
  
  twisting_probability_psitheta <- function(theta, x, ck, pkZ){
    
    BigM <- 100
    Innerprod_theta_ck <- theta*ck 
    Innerprod_theta_ck[Innerprod_theta_ck > BigM] <- BigM
    
    buff_pek <- pkZ*exp(Innerprod_theta_ck)
    
    out <- buff_pek/(1-pkZ+buff_pek)
    
    return(out)
  }
}