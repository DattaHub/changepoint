### Data 

N <- 497
theta <- numeric(N)
theta[1:138] <- -0.18
theta[139:225] <- 0.08
theta[226:242] <- 1.07
theta[243:299] <- -0.53
theta[300:308] <- 0.16
theta[309:332] <- -0.69
theta[333:n] <- -0.16

sig2 = 0.04
Y <- theta + sqrt(sig2) * rnorm(N)
plot(Y, cex=0.5, col="gray")
lines(theta)

## Stan Implementation 

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

changepoint_code <- '
data{
  int<lower=1> N;                // number of observations
  int<lower=1> L;                // number of changepoints
  vector[N] Y;                 // data matrix of size N
  real<lower=0> V;            // variance parameter 
}
transformed data {
  real I[N];
  real P[N];
  real<lower=0,upper=1> p_n;
  real<lower=0,upper=1> u;
  p_n = L*1.0/N;
  for( j in 1:N){
    u = uniform_rng(0,1);
    I[j] = (u<p_n? 0 : 1);
  }
  P <- cumulative_sum(I);
}
parameters {    
  vector[N] theta;
  real<lower=0> sigma;
  real<lower=0> tau;
  real mu_0;
  vector[L] mu;
}

model{
  sigma ~ cauchy(0, 2.5); // the hyperpriors 
  tau ~ cauchy(0,2.5); 
  mu_0 ~ normal(0, V);
  
 for (l in 1:L){
    mu[l] ~ normal(mu_0, tau);
  }
  //The likelihood
  for(j in 1:N){
    for(l in 1:L){
      if(P[j] == l-1)
        theta[j] ~ normal(mu[l],tau); 
      }
  }
 Y ~ normal(theta,sigma);
}'


data_list <- list(N = N,
                  L = length(table(theta)),
                  Y = Y, 
                  V = 10)


cp.stan.fit = stan_model(model_code= changepoint_code, model_name="BayesCP")

# mod <- stan(model_code = changepoint_code, data = data_list, chains = 1, iter = 100)

library("parallel")
Nchains <- 1
Niter <- 2000
t_start <- proc.time()[3]
sflist <- mclapply(1:Nchains, mc.cores = Nchains, 
                   function(i) stan(fit = cp.stan.fit, 
                                    data = data_list, 
                                    pars= c("theta","sigma","tau","mu_0","mu"),
                                    seed = 42,
                                    iter=Niter,
                                    #diagnostic_file = paste("diagfile",i,".csv",sep = ""),                                    sample_file = paste("sampfile",i,".csv",sep = ""),
                                    chains = 1, chain_id = i, 
                                    refresh = -1))
t_end <- proc.time()[3]
t_elapsed <- t_end - t_start


# sflist2stanfit() : Merge a list of stanfit objects into one

cp.fit<- sflist2stanfit(sflist) 
print(cp.fit,probs = c(0.5))

list_of_draws <- extract(cp.fit)
print(names(list_of_draws))
  
 