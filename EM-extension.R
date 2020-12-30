suppressMessages(require(mixtools, quietly = T))


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
} # the function is used to pre-allocate the matrix 

# the function is used to return a vecor where each entry is the sum of each element
# at same column position but different row. The new vector is used as normalizing factor
# to compute each d{kth}
sum.row <- function(A){
  res_vec <- rep(0, ncol(A))
  for (i in 1:nrow(A)){
    res_vec <- res_vec + A[i,]  
  }
  return(res_vec)
} 



handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T, k)
{
  
  like <- 0
  for (i in 1:k){
    cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3), rgb(0,0,1,.3))
    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
  }
  
  deviance <- -2*sum(log(like))
  res      <- matrix(NA,n_iter + 1, 2+3*k)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  d_tot <- rep.row(rep(NA,length(y)),k)
  r_tot <- rep.row(rep(NA,length(y)),k)
  for (iter in 1:n_iter) {
    
    # E step
    for (j in 1:k) {d_tot[j,] <- p[j]*dnorm(y, mu[j], sigma[j])}
    for (j in 1:k) {r_tot[j,] <- d_tot[j,]/sum.row(d_tot)}
    
    # M step
    for (j in 1:k){
      r = r_tot[j,]
      p[j]     <- mean(r)
      mu[j]    <- sum(r*y)/sum(r)
      sigma[j] <-sqrt( sum(r*(y^2))/sum(r) - (mu[j])^2 )
    }
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- 0
    for (i in 1:k){like <- like + p[i]*dnorm(y, mu[i], sigma[i])} 
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag && iter==n_iter){
      hist(y, prob = T, breaks=100, 
           col = "pink", border = "white",
           main = "Bart Simpson Density", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))

      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6 )
      
     likefunction <- function(y){
       like <- 0
       for (i in 1:k){
         like <- like + p[i]*dnorm(y, mu[i], sigma[i])
       }
       return(like)
     }
      
     # add curves, colors and legend
      cols = rep(NA,k+1) 
      cols[k+1] = "orange"
      # #iterate to plot every distribution
      for (i in 1:k) {
        cols[i] = rgb(runif(1),runif(1),runif(1),.8)
        curve(p[i]*dnorm(x,mu[i],sigma[i]),
                            lwd = 3,
                            lty = 5,
                            col = cols[i],
                            add = TRUE)}
      # add the join distribution
      curve(likefunction(x),
            lwd = 4, col = "orange", add = TRUE)
      # add the legend
      legendlabels=rep(NA,k+1)
      for (i in 1:k){legendlabels[i]=paste("Distribution ", i)}
      legendlabels[k+1] = "joint "
      
      legend('topright',legend=legendlabels,
             col=cols, lty=1, cex=0.6)
      grid()
      
      #Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  #names(res) <- c("iteration","p1","p2","mu1","mu2","sigma1","sigma2","deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}


n <- 3000 # Sample size
K=6 # number of ditributions
XX <- rnormmix(n,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )

hist(XX, prob = T, col = "pink",
     border = "white", breaks = 100, 
     main = "Bart Simpson", 
     xlab = "")

hem_fit <- handmade.em(XX, 
                       p      = rep(1/K,K), 
                       mu     = c(.4,.1,-.8,-.2,.8,.1),
                       sigma  = c(.1,.1,.4,.1,.6,.8), 
                       n_iter = 150,
                       plot_flag = T,
                       k=K)
round( hem_fit$parameters, 3 )
hem_fit$deviance


# AIC ---------------------------------------------------------------------


likefunction <- function(y, hem_fit){
  n_gauss = length(hem_fit$parameters)/3
  p = hem_fit$parameters[1:n_gauss]
  mu = hem_fit$parameters[(n_gauss+1):(n_gauss*2)]
  sigma = hem_fit$parameters[(n_gauss*2+1):(n_gauss*3)]
  like <- 0
  for (i in 1:n_gauss){
    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
  }
  return(log(like))
}

ll <- likefunction(XX, hem_fit)

M <- 30
k_max <- 10
results=rep(0,k_max)
for (j in 1:M){
  n <- 5000 # Sample size
  XX <- rnormmix(n,
                 lambda = c(0.5, rep(0.1,5)),
                 mu = c(0, ((0:4)/2)-1),
                 sigma = c(1, rep(0.1,5)) )
  for (i in 1:k_max){
    hem_fit <- handmade.em(XX, 
                           p      = rep(1/i,i), 
                           mu     = runif(i,min=-1.5,max=1.5),
                           sigma  = runif(i,min=0.1,max=0.4), 
                           n_iter = 500,
                           plot_flag = T,
                           k=i)
    ll <- likefunction(XX, hem_fit)
    
    AIC <- -2*sum(ll)+2*length(hem_fit$parameters)
    #print(c(AIC,length(hem_fit$parameters)/3))
    results[length(hem_fit$parameters)/3]=results[length(hem_fit$parameters)/3]+AIC
  }
  remove(XX)
  print(which.min(results))
}



AIC <- -2*mean(ll)+2*length(hem_fit$parameters)
c(AIC,length(hem_fit$parameters))
library(stats)
?AIC
AIC(ll,k=2)


