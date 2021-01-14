suppressMessages(require(mixtools, quietly = T))
suppressMessages(require(KScorrect, quietly = T))
suppressMessages(require(caret, quietly = T))
library(gt)
#library(caret)
#library("KScorrect")

gen_distr <- function(n, M){
  distr = list()
  for (j in 1:M) distr[[j]] = rnormmix(n, lambda = c(0.5, rep(0.1,5)),
                                       mu = c(0, ((0:4)/2)-1), sigma = c(1, rep(0.1,5)) )
  return(distr)
}



#Parameters in the code ------------------------------

set.seed(42)

init_num=8 #Every function that calls the EM algorithm tries to initialiaze multiple times 
#(for a fixed number of gaussians) and keeps the parameters that gave the best results
#on the training sample (which for AIC and BIC coincides with the test sample).
#Set this to 1 for a (much) faster run of the code.

n_small=8 #number of small samples from the Bart.
n_large=8 #number of large samples from the Bart.

k_max=15 #The number of gaussians used for the MoG will be from 1 to k_max.

small_n = gen_distr(300,n_small) #One can choose how many points will be in every small sample
large_n = gen_distr(50000,n_large) #One can choose how many points will be in every large sample





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
  
  d_tot <- rep.row(rep(NA,length(y)),k)
  r_tot <- rep.row(rep(NA,length(y)),k)
  
  for (iter in 1:n_iter) {
    
    # E step (get responsibilities)
    # responsibility = proportion times the Gaussian over the 2 found parameters
    for (j in 1:k) {d_tot[j,] <- pmax(p[j]*dnorm(y, mu[j], sigma[j]),1e-30)}
    # getting optimal hidden state (since it is proportionality we need to normalize) distribution
    sum_columns=sum.columns(d_tot)
    for (j in 1:k) {r_tot[j,] <- d_tot[j,]/sum_columns} 
    
    # M step
    # here we compute the p, mean and sigma for each gaussian.
    for (j in 1:k){
      r = r_tot[j,]
      sum_r=sum(r)
      p[j]     <- sum_r/length(r)
      mu[j]    <- sum(r*y)/(sum_r)   
      sigma[j] <- sqrt(max((sum(r*(y^2))/sum_r - (mu[j])^2 ),1.0e-16))   #We don't want
      #a value of sigma too small to avoid computational errors.
     }
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

temp_result=AIC_BIC(data,n_iter=200)

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
  { print(j)
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


Cross_validation5_results=cross_validation(data,5,n_iter=200)
Cross_validation10_results=cross_validation(data,10,n_iter=200)


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

Wasserstein_score_results=Wasserstein_score(data,n_iter=200)

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
write.csv2(TO_SAVE,'listHolder3.csv')
print("finished")

<<<<<<< HEAD
listHolder
typeof(length(listHolder$AIC))
######################

get_df <- function(elm_list){
  n = length(elm_list)
  res_matrix = matrix(NA, nrow = n, ncol = 4)
  for (e in 1:n){
    if (e <= n/2) res_matrix[e,] = append(elm_list[[e]], c('small_n'),after=0)
    else res_matrix[e,] = append(elm_list[[e]], c('big_n'),after=0)
  }
  df <- as.data.frame(res_matrix)
  colnames(df) <- c("sample_size","1st", "2nd", "3rd")
  return(df)
}

plt_table <- function(df,name){
  df %>% 
    gt() %>% 
    tab_header(title = name)%>%
    tab_style(style=cell_text(color='pink'), locations=cells_title(c("title")))%>%
    tab_style(style=cell_fill(color='seagreen4'), locations=cells_title(c("title"))) %>%
    
    tab_style(style=cell_fill(color='gold'), locations=cells_column_labels(columns = vars("1st"))) %>%
    data_color(columns = vars("1st"), colors = c("gold")) %>%
    
    tab_style(style=cell_fill(color='azure'), locations=cells_column_labels(columns = vars("2nd"))) %>%
    data_color(columns = vars("2nd"), colors = c("azure3")) %>%
    
    tab_style(style=cell_fill(color='sienna3'),locations=cells_column_labels(columns = vars("3rd"))) %>%
    data_color(columns = vars("3rd"), colors = c("sienna3")) %>%
    
    tab_style(style = cell_text(weight = "bold"), locations=cells_column_labels(columns = vars(sample_size)))
    
}

m_names = c("AIC","BIC","SamS_30","SamS_50","SamS_70","CV_5","CV_10","Wesserstein")
for (l in 1:length(listHolder)) print(plt_table(get_df(listHolder[[l]]),m_names[l]))

=======
>>>>>>> refs/remotes/origin/main
