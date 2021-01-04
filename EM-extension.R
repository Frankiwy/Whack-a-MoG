suppressMessages(require(mixtools, quietly = T))
set.seed(1235)

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



handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T, k, m){
  
  #function to compute the likelihood function
  likefunction <- function(y){
    like <- 0
    for (i in 1:k){
      like <- like + p[i]*dnorm(y, mu[i], sigma[i])
    }
    return(like)
  }
  
  like <- likefunction(y) # compute likelihood
  deviance <- -2*sum(log(like)) # compute deviance
  
  res      <- matrix(NA,n_iter + 1, 2+3*k)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  d_tot <- rep.row(rep(NA,length(y)),k)
  r_tot <- rep.row(rep(NA,length(y)),k)
  for (iter in 1:n_iter) {
    
    # E step (get responsibilities)
    # responsibility = proportion times the Gaussian over the 2 found parameters
    for (j in 1:k) {d_tot[j,] <- p[j]*dnorm(y, mu[j], sigma[j])}
    # getting optimal hidden state (since it is proportionality we need to normalize) distribution
    for (j in 1:k) {r_tot[j,] <- d_tot[j,]/sum.row(d_tot)} 
    
    # M step
    # here we compute the pi, mean sigma for each responsibility
    for (j in 1:k){
      r = r_tot[j,]
      p[j]     <- mean(r)
      mu[j]    <- sum(r*y)/sum(r)
      sigma[j] <-sqrt( sum(r*(y^2))/sum(r) - (mu[j])^2 ) +1.0e-8
    }
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- likefunction(y) # update likelihood 
    deviance <- -2*sum( log(like) )# update deviance
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag && iter==n_iter){
      hist(y, prob = T, breaks=100, 
           col = "pink", border = "white",
           main = paste("Bart Simpson Density, M=",m), xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))

      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6 )
      
     #likefunction <- function(y){
     #  like <- 0
     #  for (i in 1:k){
     #    like <- like + p[i]*dnorm(y, mu[i], sigma[i])
     #  }
     #  return(like)
     #}
      
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
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma),
              deviance = deviance, 
              res = res)
  return(out)
}


#n <- 3000 # Sample size
##K=6 # number of ditributions
#XX <- rnormmix(n,
#               lambda = c(0.5, rep(0.1,5)),
#               mu = c(0, ((0:4)/2)-1),
#               sigma = c(1, rep(0.1,5)) )
#
#hist(XX, prob = T, col = "pink",
#     border = "white", breaks = 100, 
#     main = "Bart Simpson", 
#     xlab = "")
#
#hem_fit <- handmade.em(XX, 
#                       p      = rep(1/K,K), 
#                       mu     = c(.4,.1,-.8,-.2,.8,.1),
#                       sigma  = c(.1,.1,.4,.1,.6,.8), 
#                       n_iter = 150,
#                       plot_flag = T,
#                       k=K)
#round( hem_fit$parameters, 3 )
#hem_fit$deviance


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

#ll <- likefunction(XX, hem_fit)

M <- 10
k_max <- 10
AIC_results=rep(0,k_max)
AIC_model= rep(0,k_max)
BIC_results=rep(0,k_max)
BIC_model= rep(0,k_max)


gen_distr <- function(n, M){
  distr = list()
  for (j in 1:M) distr[[j]] = rnormmix(n, lambda = c(0.5, rep(0.1,5)),
                                        mu = c(0, ((0:4)/2)-1), sigma = c(1, rep(0.1,5)) )
  return(distr)
}

small_n = gen_distr(750,10)
big_n = gen_distr(3000,10)
init_num = 5

set.seed(1235)
for (j in 1:M){
  XX <- small_n[[j]]
  n <- length(XX) # MANCA QUESTO
  for (i in 1:k_max){
    
    ll_sum <- -Inf
    
    for (w in 1:init_num){
      hem_fit_temp <- handmade.em(XX, 
                             p      = rep(1/i,i), 
                             mu     = runif(i,min=-1.5,max=1.5),
                             sigma  = runif(i,min=0.1,max=0.4), 
                             n_iter = 200,
                             plot_flag = F,
                             k=i,
                             m=j)
      ll_temp<- likefunction(XX, hem_fit_temp) # compute log-likelihood
      print(ll_temp)
      if (sum(ll_temp) > ll_sum) {
        ll <- ll_temp
        ll_sum <- sum(ll)
        hem_fit <- hem_fit_temp
      }
      
    }
    
    AIC <- round(-2*sum(ll)+2*length(hem_fit$parameters),2) # compute AIC
    BIC <- round(-2*sum(ll) + log(n)*length(hem_fit$parameters),2) # compute BIC
    
    AIC_results[length(hem_fit$parameters)/3]=AIC #AIC_results[length(hem_fit$parameters)/3]+AIC # compute cumulative AIC over all models
    BIC_results[length(hem_fit$parameters)/3]=BIC #BIC_results[length(hem_fit$parameters)/3]+BIC # compute cumulative BIC over all models
    
  }
  print(AIC_results)
  print(BIC_results)
  
  remove(XX)
  AIC_model[j] = which.min(AIC_results) # store best model based on AIC
  BIC_model[j] = which.min(BIC_results) # store best model based on BIC
  print(paste('AIC:',which.min(AIC_results), '& BIC:', which.min(BIC_results)))
}


# Cross-Validaton ---------------------------------------------------------

n = 1000

XX <- rnormmix(n,
               lambda = c(0.5, rep(0.1,5)),
               mu = c(0, ((0:4)/2)-1),
               sigma = c(1, rep(0.1,5)) )


cv <- function(x, k_fold){
  dataset <- x # copy the data set
  l = length(dataset) # get its length
  folds <- list() # list will contain all folds
  size <- as.integer(l/k_fold) # size of each fold
  for (n in 1:k_fold){ # iterate over all the desired folders
    fold_elements = c() # vector that will store elements
    while (length(fold_elements) < size){ # iterate until the fold size is reached
      idx_elm = sample(1:l, 1) # get index element
      fold_elements = c(fold_elements, dataset[idx_elm]) # add element to the vector
      dataset = dataset[-idx_elm] # drop element from data set
      l = length(dataset) # recompute the length
    }
    folds[[n]] = fold_elements  # add vector(representing the folder) to the folder list
  }
  #print(l)
  if (l != 0){# if there are still data points
    for (m in 1:l){ #assign them one by one to the folders until are done
      #print(m)
      folds[[m]] = c(folds[[m]], dataset[1]) # assign a data point to folder m
      dataset = dataset[-1] # delete it from the data set
    }
  }
  #print(length(data set))
  return(folds)
}

k_cv = cv(XX, 9)
#k_cv
for (i in 1:length(k_cv))print(length(k_cv[[i]]))









