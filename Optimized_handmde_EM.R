handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = F, k, m){
  
  d_tot <- rep.row(rep(NA,k),length(y))
  r_tot <- rep.row(rep(NA,k),length(y))
  
  for (iter in 1:n_iter) {
    
    # E step (get responsibilities)
    # responsibility = proportion times the Gaussian over the 2 found parameters
    for (j in 1:k) {d_tot[,j] <- pmax(p[j]*dnorm(y, mu[j], sigma[j]),1e-30)}
    # getting optimal hidden state (since it is proportionality we need to normalize) distribution
    
    r_tot=d_tot/(rowSums(d_tot,na.rm=TRUE))
    
    # M step
    # here we compute the p, mean and sigma for each gaussian.
    
    r=r_tot
    y=c(y)
    p    <- colMeans(r,na.rm=TRUE)
    col_sums=colSums(r,na.rm=TRUE) 
    mu   <- colSums(r*y)/col_sums
    sigma <- sqrt(pmax((colSums(r*(y^2),na.rm=TRUE)/col_sums- (mu)^2 ),1.0e-16))   #We don't want
    #a value of sigma too small to avoid computational errors.
    
  }
  if (plot_flag){
    
    likefunction <- function(y){
      like <- 0
      for (i in 1:k){
        like <- like + p[i]*dnorm(y, mu[i], sigma[i])
      }
      return(like)
    }
    
    hist(y, prob = T, breaks=100, 
         col = "pink", border = "white",
         main = paste("Bart Simpson Density, M=",m), xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
    
    
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
    
  }
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma))
  return(out)
  
}
