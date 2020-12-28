

handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T)
{
  # Init / 2 components only
  cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
  like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
  deviance <- -2*sum(log(like))
  res      <- matrix(NA,n_iter + 1, 8)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  
  for (iter in 1:n_iter) {
    # E step
    d1 <- p[1]*dnorm(y, mu[1], sigma[1])
    d2 <- p[2]*dnorm(y, mu[2], sigma[2])
    
    r1 <- d1/(d1 + d2)
    
    
    
    r2 <- 1 - r1
    
    # M step
    p[1]     <- mean(r1)
    mu[1]    <- sum(r1*y)/sum(r1)
    sigma[1] <-sqrt( sum(r1*(y^2))/sum(r1) - (mu[1])^2 )
    p[2]     <- 1 - p[1]
    mu[2]    <- sum((r2)*y)/sum((r2))
    sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    print(mu)
    
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
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
                       p      = c(.5,.5), 
                       mu     = c(45,55), 
                       sigma  = c(8,8), 
                       n_iter = 2,
                       plot_flag = T)
round( hem_fit$parameters, 3 )
hem_fit$deviance