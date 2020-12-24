
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
} # the function is used to pre-allocate the matrix of all the edges

sum.row <- function(A){
  res_vec <- rep(0, ncol(A))
  for (i in 1:nrow(A)){
    res_vec <- res_vec + A[i,]  
  }
  return(res_vec)
}


handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T, k)
{
  # Init / 2 components only
  
  like <- 0
  for (i in 1:k){
    cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3), rgb(0,0,1,.3))
    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
  }
  #cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
  #like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
  deviance <- -2*sum(log(like))
  res      <- matrix(NA,n_iter + 1, 2+3*k)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  d_tot <- rep.row(rep(NA,272),k)
  r_tot <- rep.row(rep(NA,272),k)
  for (iter in 1:n_iter) {
    
    # E step
    for (j in 1:k) {
      #print(c(p[j], mu[j], sigma[j]))
      d_tot[j,] <- p[j]*dnorm(y, mu[j], sigma[j])
      #print(p[j]*dnorm(y, mu[j], sigma[j]))
    }
    for (j in 1:k) {
      r_tot[j,] <- d_tot[j,]/sum.row(d_tot)
    }
    
    # M step
    for (j in 1:k){
      r = r_tot[j,]
      #print(r)
      p[j]     <- mean(r)
      mu[j]    <- sum(r*y)/sum(r)
      #print(mu[j])
      #print(sqrt( sum(r*(y^2))/sum(r) - (mu[j])^2 ))
      sigma[j] <-sqrt( sum(r*(y^2))/sum(r) - (mu[j])^2 )
      #print(sigma[j])
    }
    #p[2]     <- 1 - p[1]
    #mu[2]    <- sum((r2)*y)/sum((r2))
    #sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- 0
    for (i in 1:k){
      #cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
      like <- like + p[i]*dnorm(y, mu[i], sigma[i])
    } 
    #like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag){
      hist(y, prob = T, breaks = 30, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             col = cols[ (dnorm(y,mu[1],sigma[1]) > dnorm(y,mu[2],sigma[2])) + 1])
      
      curve(p[1]*dnorm(x,mu[1],sigma[1]) + p[2]*dnorm(x,mu[2],sigma[2]),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  names(res) <- c("iteration","p1","p2","mu1","mu2","sigma1","sigma2","deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}

data("faithful")
?faithful
hem_fit <- handmade.em(faithful$waiting, 
                       p      = c(.5,.3, .2), 
                       mu     = c(45,55,50), 
                       sigma  = c(8,8, 8), 
                       n_iter = 20,
                       plot_flag = T,
                       k=3)
round( hem_fit$parameters, 3 )
hem_fit$deviance













