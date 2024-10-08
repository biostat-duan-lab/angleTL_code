library(MASS)
library(pracma)
library(corpcor)
library(grDevices)
library(ggplot2)
library(patchwork)

simulate_ridge_risk_cov = function(Sigma, n, p, rho, cov, lambda_arr, Sigma_delta, Vt=1, Vs=1, sigma2=1){
  
  ## Monte Carlo Evaluation
  num_monte = 500
  
  a2_t = p*Vt
  a2_s = p*Vs
  gamma = p/n
  S = Sigma^(1/2)
  
  pred_risk_t = pred_risk_1D = pred_risk_2D = zeros(length(lambda_arr),1)
  for(k in 1:length(lambda_arr)){
    print(k)
    lambda = lambda_arr[k]
    pred_err_t = pred_err_1D = pred_err_2D = zeros(num_monte,1)
    for(i in 1:num_monte){
      X = pracma::randn(n,p)%*%S
      
      ###generate beta & w
      beta_w = matrix(mvrnorm(n=p, mu=c(0,0), Sigma=cov), p, 2)
      beta = matrix(beta_w[,1],p,1)
      w = matrix(beta_w[,2],p,1)
      #delta = rnorm(p, mean=0, sd = sd_delta)
      delta = mvrnorm(n=1, mu=rep(0,length(sd_delta)), Sigma = Sigma_delta)
      w = w + delta
      
      y = X%*%beta + pracma::randn(n,1)
      
      ##compute ridge
      eta = lambda * rho * sqrt(a2_t/a2_s)
      beta_t = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% t(X)%*%y
      beta_1D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*lambda*matrix(w, p, 1))
      beta_2D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*eta*matrix(w, p, 1))
      
      #inner loop, generate random test data: x_test, y_test
      x_test = pracma::randn(100,p)%*%S
      y_test = x_test%*%beta + pracma::randn(100,1)
      
      y_hat_t = x_test%*%beta_t
      y_hat_1D = x_test%*%beta_1D
      y_hat_2D = x_test%*%beta_2D
      mean((y_test - y_hat_t)^2)
      mean((y_test - y_hat_1D)^2)
      mean((y_test - y_hat_2D)^2)
      
      pred_err_t[i,1] = mean((y_test - y_hat_t)^2)
      pred_err_1D[i,1] = mean((y_test - y_hat_1D)^2)
      pred_err_2D[i,1] = mean((y_test - y_hat_2D)^2)
    }
    pred_risk_t[k] = mean(pred_err_t)
    pred_risk_1D[k] = mean(pred_err_1D)
    pred_risk_2D[k] = mean(pred_err_2D)
  }
  return(list(pred_risk_t=pred_risk_t, pred_risk_1D=pred_risk_1D, pred_risk_2D=pred_risk_2D))
}



compute_ST = function(w, t, gamma, grid_size=1e5){#grid_size=1e5
  #w,t - input spectral distribution is a mixture of point masses H = sum delta_{t_i} * w_i
  
  #the v-grid
  v = pracma::linspace(1/grid_size, 1e3, grid_size)
  z = zeros(grid_size,1)
  for(i in 1:grid_size){
    z[i] = -1/v[i] + gamma * sum(w*t/(1 + v[i]*t))
  }
  
  #find the region where z<0, this corresponds to lambda>0
  v = v[z<0]
  z = z[z<0]
  lambda = -z
  
  ind = which((lambda<10) & (lambda>1e-2))
  lambda = lambda[ind]
  v = v[ind]
  z = z[ind]
  m = v/gamma + (1/gamma-1)/z
  
  #compute m',v'
  L = length(lambda)
  v_prime = zeros(L,1)
  for(i in 1:L){
    v_prime[i] = 1 / (1/v[i]^2 - gamma*sum(w*t^2/(1 + t*v[i])^2))
  }
  m_prime = v_prime/gamma - (1/gamma-1)/z^2
  
  return(list(lambda=lambda, m=m, v=v, m_prime=m_prime, v_prime=v_prime))
}



rate = 1
plot_list = list()
order=1
sigma2 = 1
n_lambda = 100

ylabs = c("Prediction Error","Prediction Error", "Prediction Error","Prediction Error", "","", "","", "","", "","")
for(r2 in c(0.3, 0.6, 0.9)){
  for(gamma in c(0.5,2)){
    for(a in c(1,2)){
      if(a==1){
        Vt = 1^2
        Vs = 0.5^2
      }else if(a==2){
        Vt = (1/2)^2
        Vs = (0.9/2)^2
      }
      if(gamma==2){
        if(Vs==(0.9/2)^2){
          plot_list[[order]] = NA
          order = order + 1
          next
        }
        n = 50
        p = 100
        yrange = c(5,46)
        xrange = c(0.1,1)
      }else if(gamma==0.5){
        n = 50
        p = 25
        if(Vs==(0.5)^2){
          yrange = c(1.5,6.5)
        }else{
          yrange = c(1,3.1)
        }
        xrange = c(0.1,1)
      }
      
      a2_t = p*Vt
      a2_s = p*Vs
      g = pracma::linspace(1/(2*p), 1-1/(2*p), p)
      g = 1/rate*log(1/g)
      Sigma = diag(g)
      lambda = linspace(0.01, 2.5, n_lambda)^2
      t = eig(Sigma)
      w = ones(p,1)/p
      Vratio = ifelse(Vs==(0.9/2)^2, '10/9', sqrt(Vt/Vs))
      
      set.seed(123)
      sd_delta = sqrt(abs(pracma::randn(p,1)))/20
      rho_delta = 0.1
      Sigma_delta = matrix(0,length(sd_delta),length(sd_delta))
      diag(Sigma_delta) = sd_delta
      for(i in 1:nrow(Sigma_delta)){
        #print(i)
        if(i!=nrow(Sigma_delta)){
          for(j in (i+1):nrow(Sigma_delta)){
            print(j)
            Sigma_delta[i,j] = sqrt(Sigma_delta[i,i]*Sigma_delta[j,j])*rho_delta
          }
        }
      }
      Sigma_delta[lower.tri(Sigma_delta)] = Sigma_delta[upper.tri(Sigma_delta)]
      if(corpcor::is.positive.definite(Sigma_delta)==FALSE){
        Sigma_delta = corpcor::make.positive.definite(Sigma_delta)
      }
      
      set.seed(123)
      cov_raw = matrix(rbind(c(Vt, sqrt(Vt*Vs)*r2),c(sqrt(Vt*Vs)*r2, Vs)), 2, 2)
      if(corpcor::is.positive.definite(cov_raw)==FALSE){
        cov = corpcor::make.positive.definite(cov_raw)
      }else{
        cov = cov_raw
      }
      
      #empirical prediction error
      empirical_ridge = simulate_ridge_risk_cov(Sigma, n, p, r2, cov, lambda, Sigma_delta, Vt=Vt, Vs=Vs, sigma2=1)
      pred_risk_t = empirical_ridge$pred_risk_t
      pred_risk_1D = empirical_ridge$pred_risk_1D
      pred_risk_2D = empirical_ridge$pred_risk_2D
      
      #theoretical prediction error
      ST = compute_ST(w, t, gamma)
      lambda_th = ST$lambda
      m = ST$m
      v = ST$v
      m_prime = ST$m_prime
      v_prime = ST$v_prime
 
      pred_risk_th_t = (1 + (lambda_th*a2_t/gamma-1)*(1-lambda_th*v_prime/v))/(lambda_th*v)
      pred_risk_th_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s - 2*lambda_th*lambda_th*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      
      C_L = min(as.numeric(eigen(Sigma_delta)$values))*p
      C_U = max(as.numeric(eigen(Sigma_delta)$values))*p
      
      L_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s + lambda_th^2*C_L - 2*lambda_th*lambda_th*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      U_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s + lambda_th^2*C_U - 2*lambda_th*lambda_th*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      
      eta_L = (lambda_th*r2*sqrt(a2_t)*sqrt(a2_s))/(a2_s+C_L)
      eta_U = (lambda_th*r2*sqrt(a2_t)*sqrt(a2_s))/(a2_s+C_U)
      opt_eta = (eta_L+eta_U)/2
      
      L_2D = sigma2 + (lambda_th^2*a2_t + eta_L^2*a2_s + eta_L^2*C_L - 2*lambda_th*eta_L*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      U_2D = sigma2 + (lambda_th^2*a2_t + eta_U^2*a2_s + eta_U^2*C_U - 2*lambda_th*eta_U*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
 
      pred_risk_th_2D = sigma2 + (lambda_th^2*a2_t + opt_eta^2*a2_s - 2*lambda_th*opt_eta*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      
      dat = data.frame(lambda=rep(lambda,3), error=c(pred_risk_t, pred_risk_1D, pred_risk_2D), Method=c(rep('Target only empirical',length(pred_risk_t)),rep('distTL empirical',length(pred_risk_t)),rep('angleTL empirical',length(pred_risk_t))))
      dat$Method = factor(dat$Method, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
      dat_band = data.frame(lambda=rep(lambda_th,3), lower=c(pred_risk_th_t, rep(0,length(lambda_th)), L_2D), upper=c(pred_risk_th_t, rep(0,length(lambda_th)), U_2D),Color=c(rep('Target only empirical',length(pred_risk_th_t)),rep('distTL empirical',length(pred_risk_th_t)),rep('angleTL empirical',length(pred_risk_th_t))))
      dat_band$Color = factor(dat_band$Color, levels = c('Target only empirical','distTL empirical','angleTL empirical'))
      test = merge(x=dat, y=dat_band, by="lambda",all=TRUE)
      
      p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) + 
        geom_line() + 
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) + 
        xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.position = 'none') +
        xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
        ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(Vratio) ~','~gamma~'='~.(gamma)))+
        scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
        scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend")+
        scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))
      
      if(r2==0.9){
        p = ggplot(test, aes(x=sqrt(lambda),y=error, fill=Method, linetype=Method, color=Method)) +
          geom_line() +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill=Color), alpha=0.2) +
          xlim(xrange) + theme_bw() + scale_y_continuous(limits = yrange) +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14)) +
          xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
          ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(Vratio)~','~gamma~'='~.(gamma)))+
          scale_linetype_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("solid","solid","solid"))+
          scale_fill_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black",'blue','red'))+
          scale_colour_manual(label=c('Target only empirical','distTL empirical','angleTL empirical'),values=c("black", 'blue','red'), guide = "legend")
      }
      plot_list[[order]] = p
      order = order + 1
    }
  }
}

Fig4_3_by_3 = (plot_list[[3]] | plot_list[[7]] | plot_list[[11]])/
              (plot_list[[1]] | plot_list[[5]] | plot_list[[9]])/
              (plot_list[[2]] | plot_list[[6]] | plot_list[[10]])
Fig4_3_by_3

ggsave("Fig4_3_by_3.png", plot = Fig4_3_by_3, width = 16, height = 10, dpi = 300)
