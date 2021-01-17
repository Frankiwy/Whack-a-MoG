suppressMessages(require(mixtools, quietly = T))
suppressMessages(require(KScorrect, quietly = T))
suppressMessages(require(caret, quietly = T))
library(gt)
library(tidyverse)
library(readr)
library(hrbrthemes)
library(viridis)
library(forcats)
library(ggplot2)
library(gridExtra)


gen_distr <- function(n, M){
  distr = list()
  for (j in 1:M) distr[[j]] = rnormmix(n, lambda = c(0.5, rep(0.1,5)),
                                       mu = c(0, ((0:4)/2)-1), sigma = c(1, rep(0.1,5)) )
  return(distr)
}



#Parameters in the code ------------------------------

set.seed(42)

init_num=10 #Every function that calls the EM algorithm tries to initialiaze multiple times 
#(for a fixed number of gaussians) and keeps the parameters that gave the best results
#on the training sample (which for AIC and BIC coincides with the test sample).
#Set this to 1 for a (much) faster run of the code.

n_small=10 #number of small samples from the Bart.
n_large=10 #number of large samples from the Bart.

k_max=12 #The number of gaussians used for the MoG will be from 1 to k_max.

small_n = gen_distr(250,n_small) #One can choose how many points will be in every small sample
large_n = gen_distr(3000,n_large) #One can choose how many points will be in every large sample





#Useful functions and other operations------------------------------------

podium <- function(vec,size=3){
  copied_vec = vec # make copy
  result <- rep(NA,size) # initialize vec in memory
  for (i in 1:size){#iterate until the podium is full filled
    idx = which.min(copied_vec) # get the smallest
    result[i] = which.min(copied_vec)  # add it to the vec
    copied_vec[idx] = Inf # substitute it with Inf
  }
  return(result)
}



data=c(small_n,large_n) #We concatenate the lists of datasets

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
} # the function is used to pre-allocate the matrix 

# the function is used to return a vector where each entry is the sum of each element
# in that column.
sum.columns <- function(A){
  res_vec <- rep(0, ncol(A))
  for (i in 1:nrow(A)){
    res_vec <- res_vec + A[i,]  
  }
  return(res_vec)
} 



handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = F, k, m){
  
  #function to compute the likelihood function
  likefunction <- function(y){
    like <- 0
    for (i in 1:k){
      like <- like + p[i]*dnorm(y, mu[i], sigma[i])
    }
    return(like)
  }
  
  like <- likefunction(y) #compute likelihood
  like[like<1e-30]=1e-30 #sometimes we get NaN when this number is 0 or really close to it.
  #It happens rarely, most often when one of the points is really far from the distribution.
  deviance <- -2*sum(log(like)) # compute deviance
  
  res      <- matrix(NA,n_iter + 1, 2+3*k)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  d_tot <- rep.row(rep(NA,length(y)),k)
  r_tot <- rep.row(rep(NA,length(y)),k)
  
  for (iter in 1:n_iter) {
    
    # E step (get responsibilities)
    # responsibility = proportion times the Gaussian over the 2 found parameters
    for (j in 1:k) {d_tot[j,] <- pmax(p[j]*dnorm(y, mu[j], sigma[j]),1e-30)}
    # getting optimal hidden state (since it is proportionality we need to normalize) distribution
    for (j in 1:k) {r_tot[j,] <- d_tot[j,]/sum.columns(d_tot)} 
    
    # M step
    # here we compute the p, mean and sigma for each gaussian.
    for (j in 1:k){
      r = r_tot[j,]
      p[j]     <- mean(r)
      mu[j]    <- sum(r*y)/(sum(r))   
      sigma[j] <- max(sqrt( sum(r*(y^2))/sum(r) - (mu[j])^2 ),1.0e-8)   #We don't want
      #a value of sigma too small to avoid computational errors.
      if (is.nan(sigma[j])) {
        sigma[j]=0.05
        print("NaN found in sigma[j]. Value reinitialized to 0.05") 
        #Sometimes the sqrt in sigma[j] returns NaN. We decided to reinitialize that sigma
        #value to 0.05 and hope it doesn't give NaN the second time.
      }
      }
    
    # -2 x log-likelihood (a.k.a. deviance)
    like <- pmax(likefunction(y),1e-30)    # update likelihood 
    deviance <- -2*sum( log(like) )# update deviance
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag && iter==n_iter){
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
  }
  res <- data.frame(res)
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma),
              deviance = deviance, 
              res = res)
  return(out)
}



# AIC ---------------------------------------------------------------------

#This function takes in input a dataset and the results of the hem_fit function and gives
#the log_likelihood associated to the MoG with the parameters learned by the EM algorithm 
#on the given data.
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







AIC_BIC <-function(Dataset,n_iter=500){
  AIC_models=list()
  BIC_models=list()
  
  AIC_results=rep(0,k_max)
  BIC_results=rep(0,k_max)

  
  for (j in 1:length(Dataset)){ #This cycle iterates over all samples 
    XX <- Dataset[[j]]
    n <- length(XX) 
    
    for (i in 1:k_max){ #This is the cycle that tries different number of gaussians
      
      
      #We now do a for loop to find the best initialization data for a particular
      #fixed k (number of gaussians) and for a fixed dataset XX.
      
      ll_sum <- -Inf  #This is a default initialization value. Every likelihood will be better
      #than this.
      
      for (w in 1:init_num){
        
        hem_fit_temp <- handmade.em(XX, 
                               p      = rep(1/i,i), 
                               mu     = runif(i,min=-1.5,max=1.5),
                               sigma  = runif(i,min=0.1,max=0.4), 
                               n_iter = n_iter,
                               plot_flag = F,
                               k=i,
                               m=j)
        
        ll_temp<- likefunction(XX, hem_fit_temp) # compute log-likelihood on XX.
        
        if (sum(ll_temp) > ll_sum) { # if the new ll is better than the best we had before
          ll <- ll_temp    # save new ll in the ll variable
          ll_sum <- sum(ll) # also compute the sum for the next iteration
          hem_fit <- hem_fit_temp # and save all the values returned by the EM algorithm.
        }
        
      }
      n_par=length(hem_fit$parameters)-1 # if one knows n-1 values for p, the last one is forced
      #to be 1-(sum(p_i)). This means that there's one less degree of freedom in the possible
      #values of p.
      AIC <- round(-2*sum(ll)+2*n_par,2) # compute AIC
      BIC <- round(-2*sum(ll) + log(n)*n_par,2) # compute BIC
      
      AIC_results[length(hem_fit$parameters)/3]=AIC 
      BIC_results[length(hem_fit$parameters)/3]=BIC 
      
    }
    #print(AIC_results)
    #print(BIC_results)
    
    AIC_models[[j]] <- podium(AIC_results)
    BIC_models[[j]] <- podium(BIC_results)
    
    
    remove(XX)

  }
  
  out <- list(AIC = AIC_models, 
              BIC = BIC_models)
  
  return (out)
}

temp_result=AIC_BIC(data,n_iter=100)

AIC_results=temp_result$AIC
BIC_results=temp_result$BIC

# Sample splitting ---------------------------------------------------------


#sample splitting with p=0.5

sample_splitting <-function(train_size,Dataset,n_iter=500) {
  models=list()
  for (j in 1:length(Dataset))
    {
    dataset=Dataset[[j]]
    trainIndex<-createDataPartition(dataset, p = train_size, 
                        list = FALSE, 
                        times = 1)
    train <- dataset[trainIndex]
    test <- dataset[-trainIndex]
    
    
    
    results=rep(0,k_max)
    
    
    for (i in 1:k_max){
      
      ll_sum <- -Inf
      
      #THIS CYCLE KEEPS K FIXED AND IS HERE TO SELECT BETTER INITIALIZATION PARAMETERS FOR A FIXED K
      for (w in 1:init_num)
        {
        #a=runif(i,min=1,max=10)
        #a=a/sum(a)
        hem_fit_temp <- handmade.em(train, 
                                    p      = rep(1/i,i), 
                                    mu     = runif(i,min=-1.5,max=1.5),
                                    sigma  = runif(i,min=0.1,max=0.4), 
                                    n_iter = n_iter,
                                    plot_flag = F,
                                    k=i,
                                    m=j)
        ll_train<- likefunction(train, hem_fit_temp) # compute log-likelihood with TRAIN
        ll_test<- likefunction(test, hem_fit_temp) # compute log-likelihood with TEST
        if (sum(ll_train) > ll_sum) 
          {
            ll <- ll_test
            ll_sum <- sum(ll_train)
            hem_fit <- hem_fit_temp
          }
      }
      #Now ll corresponds to the log-likelihood of the test dataset with the inizialization 
      #data that had the best results on the train set
      results[length(hem_fit$parameters)/3]= -sum(ll)

    }
    
    models[[j]] <- podium(results)
   
    
  }
  return (models)
}


Sample_splitting30_results=sample_splitting(0.30,data,n_iter=100) 
Sample_splitting50_results=sample_splitting(0.50,data,n_iter=100)
Sample_splitting70_results=sample_splitting(0.70,data,n_iter=100)

# Cross-Validaton ---------------------------------------------------------
  


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


cross_validation<-function(Dataset,k_folds,n_iter){
  models=list()
  for (j in 1:length(Dataset))
  {
    dataset=Dataset[[j]]
    k_cv=cv(dataset, k_folds)
    
    results=rep(0,k_max)
    #results[length(hem_fit$parameters)/3]= -sum(ll)
    
    for (i in 1:k_max)
    { 
      
      
      k_cv_results=rep(NA,k_folds)
      
      for (fold in 1:k_folds)
        {
        test=k_cv[[fold]]
        train=do.call(c,k_cv[-fold])
        
        
        ll_sum <- -Inf
        
        #THIS CYCLE KEEPS K FIXED AND IS HERE TO SELECT BETTER INITIALIZATION PARAMETERS FOR A FIXED K
        for (w in 1:init_num)
        {
          #a=runif(i,min=1,max=10)
          #a=a/sum(a)
          hem_fit_temp <- handmade.em(train, 
                                      p      = rep(1/i,i), 
                                      mu     = runif(i,min=-1.5,max=1.5),
                                      sigma  = runif(i,min=0.1,max=0.4), 
                                      n_iter = n_iter,
                                      plot_flag = F,
                                      k=i,
                                      m=j)
          ll_train<- likefunction(train, hem_fit_temp) # compute log-likelihood with TRAIN
          ll_test<- likefunction(test, hem_fit_temp) # compute log-likelihood with TEST
          if (sum(ll_train) > ll_sum) 
          {
            ll <- ll_test
            ll_sum <- sum(ll_train)
            hem_fit <- hem_fit_temp
          }
        }
        #Now ll corresponds to the log-likelihood of the test dataset with the inizialization 
        #data that had the best results on the train set
        
        k_cv_results[fold]=sum(ll_test)
        
        
      }
      #print(k_cv_results)
      ll_mean=mean(k_cv_results)
      results[length(hem_fit$parameters)/3]= -ll_mean
      #print(ll_mean)
      #print(" ")
      
      
      
      
      
      
    }
    models[[j]] <- podium(results)
      
  }
  return (models)
}


Cross_validation5_results=cross_validation(data,5,n_iter=100)
Cross_validation10_results=cross_validation(data,10,n_iter=100)


# Wasserstein-based estimation ----------------------------------------------------


Wasserstein_score <-function(Dataset,n_iter=500) {
  models=list()
  for (j in 1:length(Dataset))
  
  {
    results=rep(0,k_max)
    
    dataset=Dataset[[j]]
    trainIndex<-createDataPartition(dataset, p = 0.5, 
                                    list = FALSE, 
                                    times = 1)
    train <- dataset[trainIndex]
    test <- dataset[-trainIndex]
    
    Quantile_to_integrate_te <-function(z) quantile(dataset,prob=z,names=FALSE)
    
    for (i in 1:k_max){
      
      ll_sum <- -Inf
      
      #THIS CYCLE KEEPS K FIXED AND IS HERE TO SELECT BETTER INITIALIZATION PARAMETERS FOR A FIXED K
      for (w in 1:init_num)
      { hem_fit_temp <- handmade.em(train, 
                                    p      = rep(1/i,i), 
                                    mu     = runif(i,min=-1.5,max=1.5),
                                    sigma  = runif(i,min=0.1,max=0.4), 
                                    n_iter = n_iter,
                                    plot_flag = F,
                                    k=i,
                                    m=j)
        ll_train<- likefunction(train, hem_fit_temp) # compute log-likelihood with TRAIN
        
        if (sum(ll_train) > ll_sum) 
        {
          ll <- ll_train
          ll_sum <- sum(ll_train)
          hem_fit <- hem_fit_temp
        }
      }
      #Now ll corresponds to the log-likelihood of the train dataset with the inizialization 
      #data that had the best results on the train set
      
      len=length(hem_fit$parameters)
      prob=hem_fit$parameters[1:(len/3)]
      means=hem_fit$parameters[(1+len/3):(2*len/3)]
      sigmas=hem_fit$parameters[(2*len/3+1):len]
      
      Quantile_to_integrate_k <-function(z) qmixnorm(z, mean=means, pro=prob, sd=sigmas,expand=0.10)
      
      function_to_integrate <-function(z) abs(Quantile_to_integrate_te(z)-Quantile_to_integrate_k(z))
      W_k=integrate(function_to_integrate,lower=0,upper=1,subdivisions = 10000)
      
      results[length(hem_fit$parameters)/3]=W_k$value
    }
    
    models[[j]] <- podium(results)
    
  }
  return (models)
  
}

Wasserstein_score_results=Wasserstein_score(data,n_iter=100)

listHolder=list(AIC=AIC_results,
                BIC=BIC_results,
                Sample_splitting30=Sample_splitting30_results,
                Sample_splitting50=Sample_splitting50_results,
                Sample_splitting70=Sample_splitting70_results,
                Cross_validation5=Cross_validation5_results,
                Cross_validation10=Cross_validation10_results,
                Wesserstein_score=Wasserstein_score_results
                )
TO_SAVE= do.call(rbind, listHolder)
write.csv2(TO_SAVE,'listHolder.csv')
print("finished")



# Plots and Results -------------------------------------------------------

listHolder2 <- read_delim("listHolder2.csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)
m_names = c("AIC","BIC","SamS_30","SamS_50","SamS_70","CV_5","CV_10","Wasserstein")

# PLOT TABLES FOR EACH METHOD 

m_names = c("AIC","BIC","SamS_30","SamS_50","SamS_70","CV_5","CV_10","Wasserstein")

list_res = list() # list where will be stored dfs
c=1
for (d in listHolder2[-1]){#iterate over every column
  mid_matrix = matrix(NA, nrow = 8, ncol = 4) # matrix where will be stored the results
  for (v in 1:length(d)){ # open a loop that goes from 1 to 8 
    get_value = eval(parse(text=d[v])) # get the value into num
    get_value = append(get_value, m_names[v],after=0) # append to the vector the method name
    mid_matrix[v,] = get_value # add the vector to the matrix
  }
  df = as.data.frame(mid_matrix)
  names(df) = c("Method",
                paste("1st(",c,")",sep=""),
                paste("2nd(",c,")",sep=""),
                paste("3rd(",c,")",sep=""))
  list_res[[c]] = df
  c = c+1
}


small_matrix = matrix(NA, nrow = 240, ncol = 2) #matrix storing top 3 for small n
big_matrix = matrix(NA, nrow = 240, ncol = 2) #matrix storing top 3 for big n
small_first_matrix = matrix(NA, nrow = 80, ncol = 2) #matrix storing top 1 for small n
big_first_matrix = matrix(NA, nrow = 80, ncol = 2) #matrix storing top 1 for small n
small_c = 1 # counter for top 3 small n
big_c = 1 # counter for top 3 big n
small_c_first = 1 # counter for top 1 small n 
big_c_first = 1 # counter for top 1 big n

for (n in 1:8){
  m = length(listHolder2) #length == 20
  get_model = listHolder2[n,2:m] # get a row of the listHolder file
  threshold = 1 # threshould to distringuish between results from small and big n
  for (v in get_model){ #iterate over each vector
    get_value = eval(parse(text=v)) # convert from string to number
    if (threshold<=10){ # if results belong to the small n
      small_first_matrix[small_c_first,] = append(get_value[1], m_names[n],after=0) # add to
      # the small_first matrix only the 1st result
      small_c_first = small_c_first +1
      for (e in get_value){#iterate over the the 3 results
        e = append(e, m_names[n],after=0) # add the method that returned those results
        small_matrix[small_c,] = e # add the vector to the small_matrix
        small_c = small_c+1
      }
    }
    else { # if results belong to big n
      big_first_matrix[big_c_first,] = append(get_value[1], m_names[n],after=0)# add to
      # the big_first matrix only the 1st result
      big_c_first = big_c_first +1
      for (e in get_value){#iterate over the the 3 results
        e = append(e, m_names[n],after=0)# add the method that returned those results
        big_matrix[big_c,] = e# add the vector to the big_matrix
        big_c = big_c+1
      }
    }
    threshold = threshold +1
  }
}  

# function to convert a matrix into a df
sort_results <- function(matrix_to_convert){
  df= as.data.frame(matrix_to_convert) #convert into df
  names(df) = c("Method", "Result") 
  
  #add sort levels
  df$Result <- factor(df$Result,levels = c("1", "2", "3", "4",
                                           "5", "6", "7", "8","9","10","11", "12"))
  return(df)
}

# function to plot histogram either for top 3 or top 1 results:
single_hist <- function(matrix_to_convert,sample_s){
  df = sort_results(matrix_to_convert) # call the sort function
  p <- df %>%
    ggplot( aes(x=Result, color=Method, fill=Method)) +
    geom_histogram(alpha=0.6, binwidth = 10, stat = 'count') +
    xlab(paste("SamSize n = ",sample_s,sep="")) +
    theme(
      axis.title.x = element_text(colour='red',face='bold'),
      axis.title.y = element_text(colour='red',face='bold'),
    ) 
  
  return(p)
}

# function to plot multiple histogram only of top 3 results:
multiple_hist <- function(matrix_to_convert,sample_s){
  df = sort_results(matrix_to_convert)
  p<- df %>%
    ggplot( aes(x=Result, color=Method, fill=Method)) +
    geom_histogram(alpha=0.6, binwidth = 5, stat = 'count') +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(1, "lines"),
      strip.text.x = element_text(size = 10,face='bold'),
      axis.title.x = element_text(hjust=.5, vjust=-1.5,colour='red',face='bold'),
      axis.title.y = element_text(hjust=.5, vjust=1.5,colour='red',face='bold'),
      plot.title = element_text(color='purple', size=12)
    ) +
    ggtitle(paste("Metric Comporison for sample size n =  ", sample_s, sep="")) +
    xlab("K") +
    ylab("Count") +
    
    facet_wrap(~Method)
  
  return(p)
}


grid.arrange(single_hist(small_first_matrix,300),
             single_hist(big_first_matrix,5000), nrow=1, ncol=2,
             top = textGrob("TOP 1st Results",gp=gpar(fontsize=15,font="bold")))


grid.arrange(single_hist(small_matrix,300),
             single_hist(big_matrix,5000), nrow=1, ncol=2,
             top = textGrob("TOP 3 Results",gp=gpar(fontsize=15,font="bold")))


multiple_hist(small_matrix, 300)
multiple_hist(big_matrix, 5000)
