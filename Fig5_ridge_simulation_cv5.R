.libPaths("/n/home11/yhan/R/x86_64-pc-linux-gnu-library/4.1")
library(MASS)
library(corpcor)
library(glmnet)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(ggplot2)
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
library(patchwork)

mean_squared_error <- function(y_est, y_test){
  mse = mean((y_est - y_test)^2)
  return(mse)
}

train_test_split <- function(x, y, size){
  n = length(y)
  train = sample(1:n, size*n)
  test = 1:n
  test = test[!test %in% train]
  x_train = x[train,]
  y_train = y[train]
  x_test = x[test,]
  y_test = y[test]
  return(list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test))
}


simulate_coef <- function(mean=c(0,0), var, rho, p,ortho = F){
  if(ortho){
    sd_delta = sqrt((1/(rho^2)-1)*var[2])
    beta = rnorm(p,mean = 0,sd = sqrt(var[1]))
    # w = rnorm(p,mean = 0,sd = sqrt(var[2]))
    delta = rnorm(p, mean=0, sd = sd_delta)
    w = beta - delta
  }
  else{
    cov_raw = matrix(rbind(c(var[1], sqrt(var[1]*var[2])*rho),c(sqrt(var[1]*var[2])*rho, var[2])), 2, 2)
    if(corpcor::is.positive.definite(cov_raw)==FALSE){
      cov = corpcor::make.positive.definite(cov_raw)
    }else{
      cov = cov_raw
    }
    beta_w = matrix(mvrnorm(n=p, mu=mean, Sigma=cov), p, 2)
    w = matrix(beta_w[,2],p,1)
    beta =  matrix(beta_w[,1],p,1)
  }
  return(return(list(beta=beta, w=w)))
}

# simulate_data <- function(n, p, beta, error=0.5){
#   x = matrix(rnorm(n*p, 0, 1), n, p)
#   y = matrix(rnorm(n, x%*%beta, error), n, 1)
#   return(list(x=x, y=y))
# }

generate_sigma <- function(p, sd_x, rho_x){
  sigma_x = matrix(0,p,p)
  diag(sigma_x) = sd_x
  for(i in 1:nrow(sigma_x)){
    if(i!=nrow(sigma_x)){
      for(j in (i+1):nrow(sigma_x)){
        sigma_x[i,j] = sqrt(sigma_x[i,i]*sigma_x[j,j])*rho_x
      }
    }
  }
  sigma_x[lower.tri(sigma_x)] = sigma_x[upper.tri(sigma_x)]
  if(corpcor::is.positive.definite(sigma_x)==FALSE){
    sigma_x = corpcor::make.positive.definite(sigma_x)
  }
  return(sigma_x)
}


simulate_data <- function(n, p, beta, error=0.5){
  sigma = generate_sigma(p, 1, 0.1)
  x = matrix(rnorm(n*p, 0, sigma), n, p)
  y = matrix(rnorm(n, x%*%beta, error), n, 1)
  # y_center = (y - mean(y))/sd(y)
  return(list(x=x, y=y))
}

simulate_data_source <- function(n, p, beta, error=0.5){
  sigma = generate_sigma(p, 1, 0.5)
  x = matrix(rnorm(n*p, 0, sigma), n, p)
  y = matrix(rnorm(n, x%*%beta, error), n, 1)
  # y_center = (y - mean(y))/sd(y)
  return(list(x=x, y=y))
}

ridge_closed_form_result <- function(x, y, lam, w, eta, x_val, y_val, standard){
  n = nrow(x)
  p = ncol(x)
  I = diag(p)
  
  
  if(n>p){
    # print(dim(t(x)%*%x + n * lam * I))
    
    if(standard==TRUE){
      beta_est = solve(t(x)%*%x + n * lam * I, tol=1e-1000) %*% t(x) %*% y
    }else{
      w_vector = matrix(w, p, 1)
      beta_est = solve(t(x)%*%x + n * lam * I, tol=1e-1000) %*% (t(x) %*% y + n * eta * w_vector)
    }
  }else{
    
    x_svd = svd(x, nu=n, nv=p)
    u = x_svd$u
    s = cbind(diag(x_svd$d), matrix(0, n, p-n))
    v = x_svd$v
    R = u %*% s
    RtR = t(R) %*% R
    RtY = t(R) %*% y
    
    if(corpcor::is.positive.definite(RtR + n * lam * I)==FALSE){
      cov = corpcor::make.positive.definite(RtR + n * lam * I)
    }else{
      cov = RtR + n * lam * I
    }
    if(standard==TRUE){
      beta_est = v %*% solve(cov, tol=1e-1000) %*% RtY
    }else{
      w_vector = matrix(w, p, 1)
      beta_est = v %*% solve(cov, tol=1e-1000) %*% (RtY + n * eta * t(v) %*% w_vector)
    }
    
  }
  mse = mean_squared_error(x_val%*%beta_est, y_val)
  return(list(beta=beta_est, mse=mse))
}

CV_ridge <- function(x_train, y_train, w, rho, var){
  mse_2D_best_tmp = list()
  mse_1D_best_tmp = list()
  mse_target_best_tmp = list()
  
  ss=c(rep(1,10),rep(2,10),rep(3,10))
  ss=sample(ss)
  n = as.integer(nrow(x_train)*2/3)
  ##############################
  ###Proposed: 2D search
  ##############################
  ##define theoretical optimal lambda & eta
  opt_lam = 1/(n * (1-rho**2))
  opt_eta = rho * opt_lam * sqrt(var[1]/var[2])
  
  ##define 2D search grid
  lam_2D_grid = c(opt_lam*1000,opt_lam*100, opt_lam*50, opt_lam*20, opt_lam*15, opt_lam*10,opt_lam*5,opt_lam*4,opt_lam*3, opt_lam*2,opt_lam,opt_lam/2,opt_lam/3,opt_lam/4,opt_lam/5, opt_lam/10, opt_lam/15, opt_lam/20, opt_lam/50, opt_lam/100)
  eta_2D_grid = c(opt_eta*1000,opt_eta*100, opt_eta*50, opt_eta*20, opt_eta*15, opt_eta*10,opt_eta*5, opt_eta*4, opt_eta*3, opt_eta*2,opt_eta,opt_eta/2,opt_eta/3,opt_eta/4,opt_eta/5, opt_eta/10, opt_eta/15, opt_eta/20, opt_eta/50, opt_eta/100, 0)
  
  ##############################
  ####TL: lambda=eta, 1D search 
  ##############################
  ##define theoretical optimal lambda
  opt_lam_1D = 1/(n * (1-rho**2))
  
  ##define 1D search grid
  lam_1D_grid = c(opt_lam_1D*1000,opt_lam_1D*100, opt_lam_1D*50, opt_lam_1D*20, opt_lam_1D*15, opt_lam_1D*10,opt_lam_1D*5, opt_lam_1D*4, opt_lam_1D*3, opt_lam_1D*2,opt_lam_1D,opt_lam_1D/2,opt_lam_1D/3,opt_lam_1D/4,opt_lam_1D/5, opt_lam_1D/10, opt_lam_1D/15, opt_lam_1D/20, opt_lam_1D/50, opt_lam_1D/100, opt_lam_1D/1000)
  
  ##############################
  ####Traditional ridge, 1D search
  ##############################
  ##define theoretical optimal lambda
  opt_lam_target = 1/n
  
  ##define 1D search grid
  lam_target_grid = c(opt_lam_target*1000,opt_lam_target*100, opt_lam_target*50, opt_lam_target*20, opt_lam_target*15, opt_lam_target*10,opt_lam_target*5,opt_lam_target*4,opt_lam_target*3, opt_lam_target*2,opt_lam_target,opt_lam_target/2,opt_lam_target/3,opt_lam_target/4,opt_lam_target/5, opt_lam_target/10, opt_lam_target/15, opt_lam_target/20, opt_lam_target/50, opt_lam_target/100, opt_lam_target/1000)
  
  for (k in 1:3) {
    # print(k)
    # ss <- sample(1:3,size=30,replace=TRUE,prob=c(1/3,1/3,1/3))
    
    x = x_train[ss!=k,]
    y = y_train[ss!=k]
    x_val = x_train[ss==k,]
    y_val = y_train[ss==k]
    
    # n = nrow(x)
    # p = ncol(x)
    
    ##begin 2D grid search
    lam_2D_best = eta_2D_best = beta_2D_best = mse_2D_best = Inf
    mse_2D = matrix(NA,length(lam_2D_grid),length(eta_2D_grid))
    for(i in 1:length(lam_2D_grid)){
      l_2D = lam_2D_grid[i]
      for(j in 1:length(eta_2D_grid)){
        e_2D = eta_2D_grid[j]
        #fit_ridge_2D = ridge_closed_form_result(lam=l_2D, eta=e_2D, x=x, y=y, x_val=x_val, y_val=y_val, w=w, standard=FALSE)
        fit_ridge_2D = ridge_closed_form_result(x=x, y=y, lam=l_2D, w=w, eta=e_2D, x_val=x_val, y_val=y_val, standard=FALSE)
        mse_2D[i,j] = fit_ridge_2D$mse
      }
    }
    mse_2D_best_tmp[[k]] = mse_2D
    # print('2D Yeah!')
    
    ##begin 1D grid search
    lam_1D_best = beta_1D_best = mse_1D_best = Inf
    mse_1D = matrix(NA,length(lam_1D_grid),1)
    for(i in 1:length(lam_1D_grid)){
      l_1D = lam_1D_grid[i]
      fit_ridge_1D = ridge_closed_form_result(lam=l_1D, eta=l_1D, x=x, y=y, x_val=x_val, y_val=y_val, w=w, standard=FALSE)
      # fit_ridge_1D = ridge_closed_form_result(n, p, v, RtR, lam=l_1D, I, RtY, w_vector, eta=l_1D, x_val=x_val, y_val=y_val, standard=FALSE)
      mse_1D[i,1] = fit_ridge_1D$mse
    }
    mse_1D_best_tmp[[k]] = mse_1D
    # print('1D Yeah!')    
    
    ##begin 1D grid search (target only)
    lam_target_best = beta_target_best = mse_target_best = Inf
    mse_target = matrix(NA,length(lam_target_grid),1)
    for(i in 1:length(lam_target_grid)){
      l_target = lam_target_grid[i]
      fit_ridge = ridge_closed_form_result(lam=l_target, eta=NULL, x=x, y=y, x_val=x_val, y_val=y_val, w=NULL, standard=TRUE)
      # fit_ridge = ridge_closed_form_result(n, p, v, RtR, lam=l_target, I, RtY, w_vector=NULL, eta=NULL, x_val=x_val, y_val=y_val, standard=TRUE)
      mse_target[i,1] = fit_ridge$mse
    }
    
    mse_target_best_tmp[[k]] = mse_target
    # print('Target Yeah!')
  }
  
  n = nrow(x_train)
  p = ncol(x_train)
  
  mse_2D_best0=mse_2D_best_tmp[[1]]+mse_2D_best_tmp[[2]]+mse_2D_best_tmp[[3]]
  #+mse_2D_best_tmp[[4]]+mse_2D_best_tmp[[5]]
  lam_2D_best = lam_2D_grid[which(mse_2D_best0==min(mse_2D_best0),arr.ind=TRUE)[1]]
  eta_2D_best = eta_2D_grid[which(mse_2D_best0==min(mse_2D_best0),arr.ind=TRUE)[2]]
  
  fit_ridge_2D = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_2D_best, w=w, eta=eta_2D_best, x_val=x_train, y_val=y_train, standard=FALSE)
  beta_2D_best = fit_ridge_2D$beta
  mse_2D_best = fit_ridge_2D$mse
  
  mse_1D_best0=mse_1D_best_tmp[[1]]+mse_1D_best_tmp[[2]]+mse_1D_best_tmp[[3]]
  #+mse_1D_best_tmp[[4]]+mse_1D_best_tmp[[5]]
  lam_1D_best = lam_1D_grid[which(mse_1D_best0==min(mse_1D_best0),arr.ind=TRUE)[1]]
  
  fit_ridge_1D = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_1D_best, w=w, eta=lam_1D_best, x_val=x_train, y_val=y_train, standard=FALSE)
  beta_1D_best = fit_ridge_1D$beta
  mse_1D_best = fit_ridge_1D$mse
  
  mse_target_best0=mse_target_best_tmp[[1]]+mse_target_best_tmp[[2]]+mse_target_best_tmp[[3]]
  #+mse_target_best_tmp[[4]]+mse_target_best_tmp[[5]]
  lam_target_best = lam_target_grid[which(mse_target_best0==min(mse_target_best0),arr.ind=TRUE)[1]]
  
  fit_ridge = ridge_closed_form_result(x=x_train, y=y_train, lam=lam_target_best, w=NULL, eta=NULL, x_val=x_train, y_val=y_train,standard=TRUE)
  beta_target_best = fit_ridge$beta
  mse_target_best = fit_ridge$mse
  
  
  return(list(mse_2D_best=mse_2D_best, mse_1D_best=mse_1D_best, mse_target_best=mse_target_best, beta_2D_best=beta_2D_best, beta_1D_best=beta_1D_best, beta_target_best=beta_target_best, lam_2D_best=lam_2D_best, eta_2D_best=eta_2D_best, lam_1D_best=lam_1D_best, lam_target_best=lam_target_best))
}

simulation <- function(iters, phos, var, p, n, n_src, ortho=FALSE){
  mse_beta = mse_w = mse_proposed_list = mse_proposed_list_test = mse_TL_list = mse_TL_list_test = mse_target_only_list = mse_target_only_list_test = lam_target_list =matrix(0, length(phos), iters)
  lam_proposed_list = eta_proposed_list = lam_TL_list = matrix(0, length(phos), iters)
  
  for(i in 1:length(phos)){
    r = phos[i]
    print(r)
    
    for(j in 1:iters){
      
      sim_coef = simulate_coef(var=var, p=p, rho=r,ortho = ortho)
      beta = sim_coef$beta
      w = sim_coef$w
      # w = w*10
      while (r+0.05<cor(beta,w) || r-0.05>cor(beta,w)) {
        sim_coef = simulate_coef(var=var, p=p, rho=r,ortho = ortho)
        beta = sim_coef$beta
        w = sim_coef$w
      }
      w_unit = apply(matrix(w,p,1), 2, function(x) x/sqrt(sum(x^2)))
      
      ## estimating w_hat
      sim_data_w = simulate_data_source(5000, p, w)
      x_w = sim_data_w$x
      y_w = sim_data_w$y
      w_hat = coef(lm(y_w~x_w))[-1]
      w_hat_unit = apply(matrix(w_hat,p,1), 2, function(x) x/sqrt(sum(x^2)))
      # print(cor(beta,w))
      
      # set.seed(1)
      sim_data = simulate_data(n=30, p, beta)
      x_train = sim_data$x
      y_train = sim_data$y
      
      sim_data = simulate_data(n=200, p, beta)
      x_test = sim_data$x
      y_test = sim_data$y
      
      # xytrain_split = train_test_split(x_norm, y_norm, size=0.6)
      # x_train = xytrain_split$x_train
      # y_train = xytrain_split$y_train
      # x_test = xytrain_split$x_test
      # y_test = xytrain_split$y_test
      
      mse_beta[i,j] = mean_squared_error((x_test %*% beta), y_test)
      mse_w[i,j] = mean_squared_error((x_test %*% w_hat), y_test)
      # 
      ##run cross validation with w
      
      fit_cv = CV_ridge(x=x_train,y=y_train, w=w_hat, rho=r, var=var)
      mse_2D_best = fit_cv$mse_2D_best
      mse_1D_best = fit_cv$mse_1D_best
      mse_target_best = fit_cv$mse_target_best
      beta_2D_best = fit_cv$beta_2D_best
      beta_1D_best = fit_cv$beta_1D_best
      beta_target_best = fit_cv$beta_target_best
      lam_2D_best = fit_cv$lam_2D_best
      eta_2D_best = fit_cv$eta_2D_best
      lam_1D_best = fit_cv$lam_1D_best
      lam_target_best = fit_cv$lam_target_best

      ### Proposed
      lam_proposed_list[i,j] = lam_2D_best
      eta_proposed_list[i,j] = eta_2D_best
      mse_proposed_list[i,j] = mse_2D_best
      mse_proposed_list_test[i,j] = mean_squared_error(x_test%*%beta_2D_best,y_test)
      
      ### TL:lambda=eta
      lam_TL_list[i,j] = lam_1D_best
      mse_TL_list[i,j] = mse_1D_best
      mse_TL_list_test[i,j] = mean_squared_error(x_test%*%beta_1D_best,y_test)
      
      ### target only (duplicated for different rho i, so only save column)
      lam_target_list[i,j] = lam_target_best
      mse_target_only_list[i,j] = mse_target_best
      mse_target_only_list_test[i,j] = mean_squared_error(x_test%*%beta_target_best,y_test)
    }
    print(mean(mse_target_only_list_test[i,]))
    print(mean(mse_TL_list_test[i,]))
    print(mean(mse_proposed_list_test[i,]))
  }
  return(list(mse_beta=mse_beta, 
              mse_w=mse_w,
              # mse_w_orth=mse_w_orth, mse_w_orth_hat=mse_w_orth_hat,
              mse_proposed_list_test=mse_proposed_list_test, 
              # mse_proposed_list_test_lm=mse_proposed_list_test_lm,
              mse_TL_list_test=mse_TL_list_test, 
              # mse_TL_list_test_lm=mse_TL_list_test_lm, 
              mse_target_only_list_test=mse_target_only_list_test,
              lam_proposed_list=lam_proposed_list, 
              eta_proposed_list=eta_proposed_list, 
              lam_TL_list=lam_TL_list, lam_target_list=lam_target_list))
  
}


iters = 200
phos_w = seq(0.3,0.95,0.05)
n = 30
n_src = 5000

all_sim_results <- list()

for (i in 1:9) {
  params_all = matrix(c(25,25,25,50,50,50,100,100,100,rep(0.5,9),rep(c(0.125,0.5,2),3)),9,3)
  vars = as.vector(params_all[i,2:3])
  p = params_all[i,1]
  var_beta = params_all[i,2]
  var_w = params_all[i,3]
  ortho = F
  
  # Run the simulation
  sim_result2 = simulation(iters=iters, phos=phos_w, var=vars, p=p, n=n, n_src=n_src, ortho=ortho)
  
  # Append the result to the list with a descriptive name
  all_sim_results[[paste0("p_", p, "_beta_", var_beta, "_w_", var_w, "_ortho_", ortho)]] <- sim_result2
}


#########################
## Visualization
#########################
id = NULL
iters = 200
xlab = seq(0.3, 0.95, 0.05)
for (i in 1:length(xlab)) {
  id = c(id, rep(xlab[i], iters))
}

p_list = list()

# Parameters setup
phos_w = seq(0.3, 0.95, 0.05)
n = 50
n_src = 5000
params_all = matrix(c(25, 25, 25, 50, 50, 50, 100, 100, 100, rep(0.5, 9), rep(c(0.125, 0.5, 2), 3)), 9, 3)
ortho = FALSE

# Loop through parameter sets
for (i in 1:9) {  # Loop over all 9 parameter sets
  # Extract the parameters
  vars = params_all[i, -1]
  p = params_all[i, 1]
  var_w = params_all[i, 3]
  var_beta = params_all[i, 2]
  ratio_p = p / n
  ratio_w = sqrt(vars[1] / vars[2])
  
  # Construct the key to retrieve the simulation result from the list
  sim_key <- paste0("p_", p, "_beta_", var_beta, "_w_", var_w, "_ortho_", ortho)
  
  # Retrieve the simulation result from the list
  sim_result_not <- all_sim_results[[sim_key]]
  
  # Now you can proceed with `sim_result_not` for further analysis or processing
  # For example:
  print(paste0("Processing simulation result for p=", p, ", beta=", var_beta, ", w=", var_w))
  
  # Prepare data for plotting
  dat_target_not = as.data.frame(as.vector(t(sqrt(sim_result_not$mse_target_only_list_test))))
  colnames(dat_target_not) = 'MSE'
  dat_target_not$Correlation = id
  dat_target_not$name = 'Target'
  
  dat_2D_not = as.data.frame(as.vector(t(sqrt(sim_result_not$mse_proposed_list_test))))
  colnames(dat_2D_not) = 'MSE'
  dat_2D_not$Correlation = id
  dat_2D_not$name = 'angleTL'
  
  dat_1D_not = as.data.frame(as.vector(t(sqrt(sim_result_not$mse_TL_list_test))))
  colnames(dat_1D_not) = 'MSE'
  dat_1D_not$Correlation = id
  dat_1D_not$name = 'distTL'
  
  dat_w_not = as.data.frame(as.vector(t(sqrt(sim_result_not$mse_w))))
  colnames(dat_w_not) = 'MSE'
  dat_w_not$Correlation = id
  dat_w_not$name = 'Source'
  
  # Combine data
  dat = rbind(dat_target_not, dat_w_not, dat_1D_not, dat_2D_not)
  dat$name = factor(dat$name, levels = c('Target', 'Source', 'distTL', 'angleTL'), labels = c('Target', 'Source', 'distTL', 'angleTL'))
  
  # Aggregate data
  dat_agg = dat %>%
    group_by(Correlation, name) %>%
    summarise(MSE = mean(MSE))
  
  # Create plot
  pic = ggplot(dat_agg, aes(x = Correlation, y = MSE, fill = name, color = name)) +
    geom_line(size = 1, aes(linetype = name)) +
    theme_bw() +
    scale_linetype_manual(name = 'Method', labels = c('Target', 'Source', 'distTL', 'angleTL'), values = c('dotted', 'dashed', 'dotdash', 'solid')) +
    scale_colour_manual(name = 'Method', labels = c('Target', 'Source', 'distTL', 'angleTL'), values = c("black", 'grey', "blue", 'red')) +
    ylim(0.5, 2) +
    xlim(0.25, 1) +
    ylab('RMSE') +
    xlab(bquote(rho)) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), plot.title = element_text(size = 16), 
          legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
    ggtitle(bquote(gamma ~ '=' ~ .(ratio_p) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(ratio_w))) +
    theme(legend.position = "none")
  
  # Store plot in list
  p_list[[i]] = pic
}

for (i in c(1, 4, 7)) {
  p_list[[i]] <- p_list[[i]] + theme(legend.position = "right")
}

# Arrange plots in 3x3 grid
(p_list[[3]] / p_list[[6]] / p_list[[9]])| 
  (p_list[[2]] / p_list[[5]] / p_list[[8]]) | 
  (p_list[[1]] / p_list[[4]] / p_list[[7]])

