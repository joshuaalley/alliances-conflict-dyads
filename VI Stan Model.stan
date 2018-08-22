// Joshua Alley
// Texas A&M University
// Varying Intercepts Logit Model 



data {
  int<lower = 0> N; // number of observations
  int<lower = 0> K; // number of IVs
  int<lower = 1> C; // number of challenger states
  int<lower = 1, upper = C> chall[N]; // challenger identifier
  int<lower = 1> T; //number of target states 
  int<lower = 1, upper = T> targ[N]; // target identifier
  int<lower = 1> D; // number of dyads
  int<lower = 1, upper = D> dyad[N]; // dyad identifier
  matrix[N, K] X; // matrix of covariates
  int<lower = 0, upper = 1> y[N]; // outcome
}

parameters {
  real alpha; // global intercept
  vector[K] beta;  // parameters on IVs

  real<lower = 0> sigma_ch; // variance hyperparameter of challenger intercepts
  real<lower = 0> sigma_tg; // variance hyperparameter of target intercepts
  real<lower = 0> sigma_dy; // variance hyperparameter for dyad intercepts
  vector[C] alpha_ch_std; // better-behaved distribution for NC parameterization
  vector[T] alpha_tg_std; // better-behaved distribution for NC parameterization
  vector[D] alpha_dy_std; // better-behaved distribution for NC parameterization
}

transformed parameters{
  // implement a non-centered parameterization of the intercepts
    vector[C] alpha_ch; // VI parameters for challengers
    vector[T] alpha_tg; // VI parameters for targets
    vector[D] alpha_dy; // VI parameters for dyads
  
  alpha_ch = 0 + sigma_ch * alpha_ch_std;
  alpha_tg = 0 + sigma_tg * alpha_tg_std;
  alpha_dy = 0 + sigma_dy * alpha_dy_std;
  
}

model {
  
  vector[N] y_hat; // linear prediction 
  
  for(i in 1:N)
  y_hat[i] = alpha + alpha_ch[chall[i]] + alpha_tg[targ[i]] + alpha_dy[dyad[i]] +  X[i] * beta; 
  
  alpha ~ normal(0, 3);
  alpha_ch_std ~ normal(0, 1);
  sigma_ch ~ normal(0, 1); // half-normal
  alpha_tg_std ~ normal(0, 1);
  sigma_tg ~ normal(0, 1); // half-normal
  alpha_dy_std ~ normal(0, 1); 
  sigma_dy ~ normal(0, 1); // half-normal
  beta ~ student_t(5, 0, 2.5);
  
  y ~ bernoulli_logit(y_hat);
}
