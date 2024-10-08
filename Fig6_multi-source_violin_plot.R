library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(ggplot2)
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
library(patchwork)

path = '/Users/tiangu/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Postdoc@Harvard/Tian Gu-shared/Angle-based Ridge/JRSSB_R1/'
iters = 1000
datas = c("multi_w_p_100_r_0.40v0.45v0.50v0.55v0.60_n_100_iter_1000_1016.RData",
          "multi_w_p_100_r_0.1v0.3v0.5v0.7v0.9_n_100_iter_1000_1016_2.RData")
titles = c('(0.4, 0.45, 0.5, 0.55,0.6)',
           '(0.1, 0.3, 0.5, 0.7, 0.9)')
numbers = c('(A)', '(B)')
p_list = list()
for(i in 1:2){
  # i=1
  load(paste0(path, datas[i]))
  title = titles[i]
  number = numbers[i]
  dat_target = as.data.frame(as.vector(t(sqrt(sim_result$mse_target_only_list_test))))
  colnames(dat_target) = 'RMSE'
  dat_target$name = 'target-only'
  
  dat_1D = as.data.frame(as.vector(t(sqrt(sim_result$mse_TL_list_test))))
  colnames(dat_1D) = 'RMSE'
  dat_1D$name = 'TL 1D'
  
  dat_2D = as.data.frame(as.vector(t(sqrt(sim_result$mse_proposed_list_test))))
  colnames(dat_2D) = 'RMSE'
  dat_2D$name = 'Proposed 2D'
  
  dat_2D_lm = as.data.frame(as.vector(t(sqrt(sim_result$mse_proposed_list_test_lm))))
  colnames(dat_2D_lm) = 'RMSE'
  dat_2D_lm$name = 'Proposed 2D_lm'
  
  dat_2D_single = as.data.frame(as.vector(t(sqrt(sim_result$mse_proposed_single_list_test))))
  colnames(dat_2D_single) = 'RMSE'
  dat_2D_single$name = 'Proposed 2D single'
  
  dat = rbind(dat_target,
              dat_2D_single, dat_2D_lm,dat_2D)
  dat$group = c(rep('target-only', iters), rep('angleTL', iters), rep('angleTL-multi1', iters), rep('angleTL-multi2', iters))
  dat$group = factor(dat$group, levels = c('target-only','angleTL', 'angleTL-multi1','angleTL-multi2'))
  dat$name = factor(dat$name, levels = c('target-only','Proposed 2D single','Proposed 2D_lm','Proposed 2D'),
                    labels = c('target-only','angleTL','angleTL-multi1','angleTL-multi2'))
  
  p = ggplot(dat, aes(x=name, y=RMSE, fill=name)) +
    scale_y_continuous(limits = c(0.1,2)) +
    geom_flat_violin(alpha = .4, trim = FALSE) +
    geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
    stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, position = position_dodge(.175)) +
    ggtitle(bquote(.(number) ~ rho ~ '=' ~ .(title))) +
    scale_fill_brewer(palette = "Dark2", name = "") +
    ylab('RMSE')+
    xlab('')+
    #facet_grid(cols = vars(group), scales = "free", space = "free") +
    # ggtitle(bquote(lambda ~ '=' ~ .(ratio_p) ~ ','~  alpha[t] ~'/' ~alpha[s] ~'=' ~ .(ratio_w)))+
    theme(text = element_text(size=13),
          strip.text = element_text(size = rel(1)),
          plot.title = element_text(size=16),
          axis.text.x=element_text(angle=10, hjust=1, size=12),
          legend.position = 'none',
          panel.background = element_rect(fill = "white",colour = "black", size = 1),
          panel.grid.major.y=element_line(color="grey",linetype=1),
          panel.grid.minor.y=element_line(color="grey", linetype=2))
  p_list[[i]] = p
}

#p_list[[1]] | p_list[[2]]
#ggsave(path = path, width = 13, height = 5, filename='multi2.png', dpi=200)

####################################
### rho=c(0.1, 0.1, 0.1, 0.1, 0.1)
###################################

path = '/Users/tiangu/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Postdoc@Harvard/Tian Gu-shared/Angle-based Ridge/JRSSB_R1/'
load(paste0(path, 'multi_w_p_50_r_0.1v0.1v0.1v0.1v0.1_iter_100_N_5000_noise.RData'))
title = '(0.1, 0.1, 0.1, 0.1, 0.1)'

dat_target = as.data.frame(as.vector(t(sim_result$mse_target)))
colnames(dat_target) = 'MSE'
dat_target$name = 'Target only'

dat_1D = as.data.frame(as.vector(t(sim_result$mse_distTL)))
colnames(dat_1D) = 'MSE'
dat_1D$name = 'Distance'

dat_multi1 = as.data.frame(as.vector(t(sim_result$mse_multi1)))
colnames(dat_multi1) = 'MSE'
dat_multi1$name = 'angleTL_multi1'

dat_multi2 = as.data.frame(as.vector(t(sim_result$mse_multi2)))
colnames(dat_multi2) = 'MSE'
dat_multi2$name = 'angleTL_multi2'

dat_w = as.data.frame(as.vector(t(sim_result$mse_w_ensemble)))
colnames(dat_w) = 'MSE'
dat_w$name = 'angleTL'


dat = rbind(dat_target,
            #dat_1D, 
            dat_w,dat_multi1, dat_multi2)
dat$group = c(rep('Target', 100), 
              #rep('Distance', 100), 
              rep('angleTL', 100), rep('angleTL-multi1', 100), rep('angleTL-multi2', 100))
dat$group = factor(dat$group, levels = c('Target', 'angleTL', 'angleTL-multi1', 'angleTL-multi2'))
dat$name = factor(dat$name, levels = c('Target only', 'angleTL', 'angleTL_multi1','angleTL_multi2'),
                  labels = c('target','angleTL','angleTL-multi1','angleTL-multi2'))

p3 = ggplot(dat, aes(x=name, y=sqrt(MSE), fill=name)) +
  scale_y_continuous(limits = c(0, 2)) +
  geom_flat_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, position = position_dodge(.175)) +
  ggtitle(bquote('(C)' ~ rho ~ '=' ~ .(title))) +
  scale_fill_brewer(palette = "Dark2", name = "") +
  xlab('')+ylab('RMSE')+
  theme_bw()+
  #facet_grid(cols = vars(group), scales = "free", space = "free") +
  theme(text = element_text(size=13),
        strip.text = element_text(size = rel(1)),
        plot.title = element_text(size=16),
        axis.text.x=element_text(angle=10, hjust=1, size=12),
        legend.position = 'none',
        panel.background = element_rect(fill = "white",colour = "black", size = 1),
        panel.grid.major.y=element_line(color="grey",linetype=1),
        panel.grid.minor.y=element_line(color="grey", linetype=2))

p_list[[3]] = p3


#############################
### K=50 sources
############################
path = '/Users/tiangu/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Postdoc@Harvard/Tian Gu-shared/Angle-based Ridge/JRSSB_R1/'
load(paste0(path, 'multi_w_p400_K50.RData'))

dat_target = as.data.frame(as.vector(t(sim_result$mse_target_only_list_test)))
colnames(dat_target) = 'MSE'
dat_target$name = 'Target only'

dat_1D = as.data.frame(as.vector(t(sim_result$mse_TL_list_test)))
colnames(dat_1D) = 'MSE'
dat_1D$name = 'TL 1D'

dat_multi1 = as.data.frame(as.vector(t(sim_result$mse_angleTL_multi1)))
colnames(dat_multi1) = 'MSE'
dat_multi1$name = 'angleTL_multi1'

dat_multi2 = as.data.frame(as.vector(t(sim_result$mse_angleTL_multi2)))
colnames(dat_multi2) = 'MSE'
dat_multi2$name = 'angleTL_multi2'

dat_w = as.data.frame(as.vector(t(sim_result$mse_best_single)))
colnames(dat_w) = 'MSE'
dat_w$name = 'angleTL'


dat = rbind(dat_target,
            dat_w, dat_multi1, dat_multi2)
dat$group = c(rep('Target', 100), rep('angleTL', 100), rep('angleTL-multi1', 100), rep('angleTL-multi2', 100))
dat$group = factor(dat$group, levels = c('Target', 'angleTL', 'angleTL-multi1', 'angleTL-multi2'))
dat$name = factor(dat$name, levels = c('Target only','angleTL','angleTL_multi1','angleTL_multi2'),
                  labels = c('target','angleTL','angleTL-multi1','angleTL-multi2'))

p4 = ggplot(dat, aes(x=name, y=sqrt(MSE), fill=name)) +
  scale_y_continuous(limits = c(0.25, 1.2)) +
  geom_flat_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, position = position_dodge(.175)) +
  ggtitle('50 source estimates') +
  ggtitle(bquote('(D) K=50 with ' ~ rho ~'between (0.1, 0.9)')) +
  scale_fill_brewer(palette = "Dark2", name = "") +
  xlab('')+ylab('RMSE')+
  theme_bw()+
  #facet_grid(cols = vars(group), scales = "free", space = "free") +
  theme(text = element_text(size=13),
        strip.text = element_text(size = rel(1)),
        plot.title = element_text(size=16),
        axis.text.x=element_text(angle=10, hjust=1, size=12),
        legend.position = 'none',
        panel.background = element_rect(fill = "white",colour = "black", size = 1),
        panel.grid.major.y=element_line(color="grey",linetype=1),
        panel.grid.minor.y=element_line(color="grey", linetype=2))


p_list[[4]] = p4


#p_list[[3]] | p_list[[4]]
(p_list[[1]] | p_list[[2]]) / (p_list[[3]] | p_list[[4]] )









