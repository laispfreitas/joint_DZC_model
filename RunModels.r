# This code runs the models introduced in the paper
# Simultaneous {\em Aedes}-borne diseases outbreaks: a joint model for the total cases and disease-specific probability of dengue, Zika, and chikungunya
# by Alexandra M. Schmidt, La??s P. Freitas, Oswaldo G. Cruz and Marilia S. Carvalho
# 22/02/2022

rm(list=ls())


library(rstan)
require(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#loading a set of artificial data generated from model M4, the covariates are randomly sampled from normal distributions

load("SampleExample.RData")

#storing the variables that will be used in the models

total<-sample$totalS
y<-sample$disrepT
E<-sample$E
x<-sample$x
N<-sample$N
N_edges<-sample$n_edges
node1<-sample$node1
node2<-sample$node2

###############################################
#
# modelos Poisson-Multinomial M0
#
##############################################
mod0 <- stan("dzc_mod0_dens.stan", 
             data=list(N=N,p=3,
                       N_edges=N_edges,node1=node1,
                       x=x,K=3,
                       node2=node2,total=total,
                       y=y,E=E), 
             chains=3, iter=10000,warmup=3000,thin=14,
             control = list(adapt_delta = 0.8,max_treedepth=15),init="random",init_r=0.50,
             pars = c("mu","alpha","beta","phi","sigma_phi","pii","log_lik","yfit","totfit"),include=TRUE)

saveRDS(mod0, "dzc_mod0_dens.rds")

traceplot(mod0,pars=c("alpha","beta"))

#####################################
#
# Model 1
#
#####################################

mod1 <- stan("dzc_mod1_dens.stan", 
             data=list(N=N,p=3,
                       N_edges=N_edges,node1=node1,
                       x=x,K=3,
                       node2=node2,total=total,
                       y=y,E=E), 
             chains=3, iter=10000,warmup=3000,thin=14,
             control = list(adapt_delta = 0.8,max_treedepth=15),init="random",init_r=0.20,
             pars = c("mu","alpha","beta","muS","gamma","sigma_phi","pii","log_lik","yfit","totfit"),include=TRUE)

saveRDS(mod1, "dzc_mod1_dens.rds")

traceplot(mod1,pars=c("alpha","beta"))

#####################################
#
# Model 2
#
#####################################


mod2 <- stan("dzc_mod2_dens.stan", 
             data=list(N=N,p=3,
                       N_edges=N_edges,node1=node1,
                       x=x,K=3,
                       node2=node2,total=total,
                       y=y,E=E), 
             chains=3, iter=10000,warmup=3000,thin=14,
             control = list(adapt_delta = 0.8,max_treedepth=15),init="random",init_r=0.20,
             pars = c("mu","alpha","gamma", "beta","muS","U","sigma_phi","sigmau","pii","log_lik","yfit","totfit"),include=TRUE)

traceplot(mod2,pars=c("alpha","beta"))

#####################################
#
# Model 3
#
#####################################


mod3 <- stan("dzc_mod3_dens.stan", 
             data=list(N=N,p=3,
                       N_edges=N_edges,node1=node1,
                       x=x,K=3,
                       node2=node2,total=total,
                       y=y,E=E), 
             chains=3, iter=10000,warmup=3000,thin=14,
             control = list(adapt_delta = 0.8,max_treedepth=15),init="random",init_r=0.20,
             pars = c("mu","beta","alpha","muS","sigma_phi","theta","pii","tau_theta","log_lik","yfit","totfit"),include=TRUE)

traceplot(mod3,pars=c("alpha","beta"))

#####################################
#
# Model 4
# Poisson-Multinomial Wishart
#
#####################################


mod4 <- stan("dzc_mod4_densWishart.stan", 
             data=list(N=N,p=3,
                       N_edges=N_edges,node1=node1,zeros=rep(0,3),SS=diag(1,3),
                       x=x,K=3,
                       node2=node2,total=total,
                       y=y,E=E), 
             chains=3, iter=10000,warmup=3000,thin=14,
             control = list(adapt_delta = 0.8,max_treedepth=15),init="random",init_r=0.70,
             pars = c("mu","beta","alpha","muS","sigma_phi","theta","pii","Sigma","log_lik","yfit","totfit"),include=TRUE)


traceplot(mod4,pars=c("alpha","beta"))
#####################################
#
# Model 5
# Poisson-Multinomial MCAR Wishart Cholesky
#
#####################################

mod5 <- stan("dzc_mod5_MCAR_densWishart.stan", 
                data=list(N=N,p=3,
                          N_edges=N_edges,node1=node1,ones=rep(1,3),SS=diag(1,3),
                          x=x,K=3,
                          node2=node2,total=total,
                          y=y,E=E), 
                chains=3, iter=40000,warmup=20000,thin=20,
                control = list(adapt_delta = 0.95,max_treedepth=15),init="random",init_r=0.70,
                pars = c("mu","beta","alpha","muS","pii","Sigma","R","log_lik","yfit","totfit"),include=TRUE)

traceplot(mod5,pars=c("alpha","beta"))
