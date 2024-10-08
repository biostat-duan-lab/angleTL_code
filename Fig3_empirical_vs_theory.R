library(MASS)
library(pracma)
library(corpcor)
library(ggplot2)
library(patchwork)

simulate_ridge_risk_cov = function(Sigma, n, p, rho, lambda_arr, Vt=1, Vs=1, sigma2=1){
  
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

      cov_raw = matrix(rbind(c(Vt, sqrt(Vt*Vs)*rho),c(sqrt(Vt*Vs)*rho, Vs)), 2, 2)
      if(corpcor::is.positive.definite(cov_raw)==FALSE){
        cov = corpcor::make.positive.definite(cov_raw)
      }else{
        cov = cov_raw
      }
      beta_w = matrix(mvrnorm(n=p, mu=c(0,0), Sigma=cov), p, 2)
      beta = matrix(beta_w[,1],p,1)
      w = matrix(beta_w[,2],p,1)
      
      y = X%*%beta + pracma::randn(n,1)
      
      ##compute ridge
      eta = lambda * rho * sqrt(a2_t/a2_s)
      if(n>p){
        beta_t = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% t(X)%*%y
        beta_1D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*lambda*matrix(w, p, 1))
        beta_2D = solve(t(X)%*%X +n*lambda*eye(p), tol=1e-100) %*% (t(X)%*%y + n*eta*matrix(w, p, 1))
      }else{
        x_svd = svd(X, nu=n, nv=p)
        u = x_svd$u
        s = cbind(diag(x_svd$d), matrix(0, n, p-n))
        v = x_svd$v
        w_vector = matrix(w, p, 1)
        R = u %*% s
        RtR = t(R) %*% R
        RtY = t(R) %*% y
        I = diag(p)
        x_svd = svd(X, nu=n, nv=p)
        u = x_svd$u
        s = cbind(diag(x_svd$d), matrix(0, n, p-n))
        v = x_svd$v
        w_vector = matrix(w, p, 1)
        R = u %*% s
        RtR = t(R) %*% R
        RtY = t(R) %*% y
        I = diag(p)
        
        beta_t = v %*% solve(RtR + n * lambda * I, tol=1e-100) %*% RtY
        beta_1D = v %*% solve(RtR + n * lambda * I, tol=1e-100) %*% (RtY + n * lambda * t(v) %*% w_vector)
        beta_2D = v %*% solve(RtR + n * lambda * I, tol=1e-100) %*% (RtY + n * eta * t(v) %*% w_vector)
      }
      
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



compute_ST = function(w, t, gamma, grid_size=1e5){
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



###### plot
rate = 1
sigma2 = 1
plot_list = list()
order=1
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
        a2_t = p*Vt
        a2_s = p*Vs
        Vratio = Vs/Vt
        yrange = c(5,30)
        xrange = c(0.1,1)
      }else if(gamma==0.5){
        n = 50
        p = 25
        a2_t = p*Vt
        a2_s = p*Vs
        Vratio = Vs/Vt
        if(Vs==(0.5)^2){
          yrange = c(1.5,6.5)
        }else{
          yrange = c(1,3.1)
        }
        xrange = c(0.1,1)
      }

      g = pracma::linspace(1/(2*p), 1-1/(2*p), p)
      g = 1/rate*log(1/g)
      Sigma = diag(g)
      lambda = linspace(0.001, 1.6, n_lambda)^2
      t = eig(Sigma)
      w = ones(p,1)/p
      
      #empirical prediction error
      empirical_ridge = simulate_ridge_risk_cov(Sigma, n, p, r2, lambda, Vt=Vt, Vs=Vs, sigma2=1)
      pred_risk_t = empirical_ridge$pred_risk_t
      pred_risk_1D = empirical_ridge$pred_risk_1D
      pred_risk_2D = empirical_ridge$pred_risk_2D
      
      #theoretical prediction error
      ST = compute_ST(w, t, gamma)
      lambda_th = ST$lambda
      opt_lam = sigma2 * gamma/(a2_t * (1-r2**2))
      m = ST$m
      v = ST$v
      m_prime = ST$m_prime
      v_prime = ST$v_prime
      eta = r2 * lambda_th * sqrt(a2_t/a2_s)
      
      pred_risk_th_t = (1 + (lambda_th*a2_t/gamma-1)*(1-lambda_th*v_prime/v))/(lambda_th*v)
      pred_risk_th_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s - 2*lambda_th*lambda_th*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      pred_risk_th_2D = sigma2 + (lambda_th^2*a2_t + eta^2*a2_s - 2*lambda_th*eta*r2*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
      
      dat = data.frame(lambda=c(rep(lambda,3), rep(lambda_th,3)),
                       error=c(pred_risk_t, pred_risk_1D, pred_risk_2D, pred_risk_th_t, pred_risk_th_1D, pred_risk_th_2D),
                       Color=c(rep('target-only empirical',length(pred_risk_t)),rep('distTL empirical',length(pred_risk_t)),rep('angleTL empirical',length(pred_risk_t)),
                               rep('target-only theory',length(pred_risk_th_t)),rep('distTL theory',length(pred_risk_th_t)),rep('angleTL theory',length(pred_risk_th_t))))
      dat$Color = factor(dat$Color, levels = c('target-only empirical','target-only theory','distTL empirical','distTL theory','angleTL empirical','angleTL theory'))
      
      Vratio = ifelse(Vs==(0.9/2)^2, '10/9', sqrt(Vt/Vs))
      ggp = ggplot(dat, aes(x=sqrt(lambda), y=error, linetype=Color, color=Color)) +
        geom_line(size=1) +
        theme_bw() +
        scale_y_continuous(limits = yrange) +
        scale_x_continuous(limits = xrange) +
        geom_vline(xintercept = sqrt(opt_lam)) +
        xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14),
              legend.position = "none") +
        ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(Vratio)~','~gamma~'='~.(gamma)))+
        scale_linetype_manual(name='Method', labels=c('target-only empirical','target-only theory','distTL empirical','distTL theory','angleTL empirical','angleTL theory'), 
                              values=c("solid", "dotdash","solid", "dotdash","solid", "dotdash"))+
        scale_colour_manual(name='Method', labels=c('target-only empirical','target-only theory','distTL empirical','distTL theory','angleTL empirical','angleTL theory'), 
                            values=c("grey","black", "skyblue",'blue','pink','red'))
      
      if(r2==.9){
        ggp = ggplot(dat, aes(x=sqrt(lambda), y=error, linetype=Color, color=Color)) +
          geom_line(size=1) +
          theme_bw() +
          scale_y_continuous(limits = yrange) +
          scale_x_continuous(limits = xrange) +
          geom_vline(xintercept = sqrt(opt_lam)) +
          xlab(expression(sqrt(lambda))) + ylab(ylabs[order]) +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14)) +
          ggtitle(bquote(rho ~ '=' ~ .(r2) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(Vratio)~','~gamma~'='~.(gamma)))+
          scale_linetype_manual(name='Method', labels=c('Target only empirical','Target only theory','distTL empirical','distTL theory','angleTL empirical','angleTL theory'), 
                                values=c("solid", "dotdash","solid", "dotdash","solid", "dotdash"))+
          scale_colour_manual(name='Method', labels=c('Target only empirical','Target only theory','distTL empirical','distTL theory','angleTL empirical','angleTL theory'), 
                              values=c("grey","black", "skyblue",'blue','pink','red'))
      }
      plot_list[[order]] = ggp
      order = order + 1
    }
  }
}
Fig3_3_by_3 = (plot_list[[3]] | plot_list[[7]] | plot_list[[11]])/
              (plot_list[[1]] | plot_list[[5]] | plot_list[[9]])/
              (plot_list[[2]] | plot_list[[6]] | plot_list[[10]])
Fig3_3_by_3

ggsave("Fig3_3_by_3.png", plot = Fig3_3_by_3, width = 16, height = 10, dpi = 300)

