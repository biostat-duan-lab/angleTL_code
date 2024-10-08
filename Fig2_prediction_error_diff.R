library(pracma)

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

rho = 0.9
rate = 1
n_lambda = 100
sigma2 = 1
Vt = 1^2
Vs = 0.5^2
Vratio = Vs/Vt
n=50
dat_diff_all = gamma_all = NULL
for(p in 5:500){
  print(p)
  a2_t = p*Vt
  a2_s = p*Vs
  g = pracma::linspace(1/(2*p), 1-1/(2*p), p)
  g = 1/rate*log(1/g)
  Sigma = diag(g)
  lambda = linspace(0.001, 1.6, n_lambda)^2
  t = eig(Sigma)
  w = ones(p,1)/p
  gamma = p/n
  
  #theoretical prediction error
  ST = compute_ST(w, t, gamma)
  lambda_th = ST$lambda
  m = ST$m
  v = ST$v
  m_prime = ST$m_prime
  v_prime = ST$v_prime
  eta = rho * lambda_th * sqrt(a2_t/a2_s)
  library(ggplot2)
  pred_risk_th_t = (1 + (lambda_th*a2_t/gamma-1)*(1-lambda_th*v_prime/v))/(lambda_th*v)
  pred_risk_th_1D = sigma2 + (lambda_th^2*a2_t + lambda_th^2*a2_s - 2*lambda_th*lambda_th*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
  pred_risk_th_2D = sigma2 + (lambda_th^2*a2_t + eta^2*a2_s - 2*lambda_th*eta*rho*sqrt(a2_t*a2_s) - lambda_th*gamma*sigma2)*(v-lambda_th*v_prime)/(gamma*(lambda_th*v)^2) + sigma2*(1/lambda_th/v-1)
  diff_1D_2D = pred_risk_th_1D - pred_risk_th_2D
  
  opt_lam_2D = sigma2 * gamma/(a2_t * (1-rho**2))
  opt_lam_1D = sigma2 * gamma/(a2_t + a2_s - 2*rho*sqrt(a2_t*a2_s))
  dat = data.frame(lambda=rep(lambda_th,2),
                   error=c(pred_risk_th_1D, pred_risk_th_2D),
                   Color=c(rep('distTL theory',length(pred_risk_th_t)),rep('angleTL theory',length(pred_risk_th_t))))
  dat$Color = factor(dat$Color, levels = c('distTL theory','angleTL theory'))
  
  lo_1D <- loess(error~sqrt(lambda), data=dat[dat$Color=='distTL theory',])
  lo_2D <- loess(error~sqrt(lambda), data=dat[dat$Color=='angleTL theory',])
  dat_diff = predict(lo_1D, data.frame(lambda=sqrt(opt_lam_1D))) - predict(lo_2D, data.frame(lambda=sqrt(opt_lam_2D)))
  dat_diff_all = c(dat_diff_all, dat_diff)
  gamma_all = c(gamma_all, gamma)
}

dat_diff = data.frame(cbind(gamma_all, dat_diff_all))
colnames(dat_diff) = c('gamma','diff')
ggplot(dat_diff, aes(x=gamma_all, y=dat_diff_all)) +
  geom_smooth(method = "loess", span=0.08, level=0.90, method.args=list(degree =1, family = 'gaussian'), color='black')+
  #geom_line(size=1) + 
  theme_bw() +
  xlab(expression(gamma)) + ylab('Prediction Error Difference') +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),plot.title = element_text(size=20),legend.title = element_text(size=14),legend.text = element_text(size=14))+
  ggtitle(bquote(rho ~ '=0.9' ~ ','~  alpha[t] ~ '/' ~ alpha[s] ~ '=2'))


