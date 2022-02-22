data {
  int<lower=0> p;
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=0> total[N];
  int<lower=0> y[N,3]; // count outcomes of each of the diseases
  vector<lower=0> [N] E;           // expected number of cases
  int<lower=1> K;                 // num covariates
  matrix[N, K] x;// design matrix
  vector [p] zeros;
  cov_matrix[p] SS;
}

  


parameters {
  
  vector[(K+1)] beta;
  matrix[(K+1),(p-1)] alpha;       // covariates probability
  real<lower=0> tau_phi;        // precision spatial effect
  vector [N] phi;         // spatial effects
  matrix [N,3] theta;         // independent effects
  cov_matrix[3] Sigma;
  
  
}


transformed parameters{
  
  matrix<lower=0> [N,2] lambda;
  matrix<lower=0,upper=1> [N,3] pii;
  real<lower=0> sigma_phi = inv(sqrt(tau_phi)); // convert precision to sigma
  vector [N] muS;
  vector<lower=0> [N] mu;
  
  
  muS=sigma_phi*phi;
  for(i in 1:N){
    mu[i]=exp(beta[1]+x[i,1]*beta[2]+x[i,2]*beta[3]+x[i,3]*beta[4]+muS[i]+theta[i,1]);
    for(c in 1:2){
      lambda[i,c]=exp(alpha[1,c]+x[i,1]*alpha[2,c]+x[i,2]*alpha[3,c]+x[i,3]*alpha[4,c]+muS[i]+theta[i,(c+1)]);
    }
    pii[i,1]=1/(1+sum(lambda[i,1:2]));
    pii[i,2]=lambda[i,1]/(1+sum(lambda[i,1:2]));
    pii[i,3]=lambda[i,2]/(1+sum(lambda[i,1:2]));
  }
  
}


model {
  
  // joint likelihood function
  for(i in 1:N){
      total[i] ~ poisson((E[i]*mu[i]));
      y[i,1:3] ~ multinomial(pii[i,1:3]');
  }
  
  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  
  for(i in 1:N){
    theta[i,]~multi_normal(zeros',Sigma);
  }
  
  Sigma ~ inv_wishart(10,SS);
  
  
  for(i in 1:(K+1)){
   beta[i]~normal(0,2);
   for(j in 1:(p-1))
    alpha[i,j] ~ normal(0,2);
  }  
    
    
  tau_phi ~ cauchy(0, 1);
  
  
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)
  
}

generated quantities {
  
  vector [N] log_lik;
  int<lower=0> totfit [N];
  int<lower=0> yfit[N,3];
  
  for(i in 1:N){
    log_lik[i]=poisson_lpmf(total[i] | E[i]*mu[i]) + multinomial_lpmf(y[i,1:3]|pii[i,1:3]');
    totfit[i]=poisson_rng(E[i]*mu[i]);
    yfit[i,]=multinomial_rng(pii[i,1:3]',total[i]);
    
  }


}
