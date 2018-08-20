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
  matrix[N, K] X; // matrix of covariates
  int<lower = 0, upper = 1> y[N]; // outcome
}

parameters {
  real alpha; // global intercept
  vector[K] beta;  // parameters on IVs
  vector[C] alpha_ch; // VI parameters for challengers
  real<lower = 0> sigma_ch; // variance hyperparameter of challenger intercepts
  vector[T] alpha_tg; // VI parameters for targets
  real<lower = 0> sigma_tg; // variance hyperparameter of target intercepts
}

model {
  
  vector[N] y_hat; // linear prediction 
  
  for(i in 1:N)
  y_hat[i] = alpha + alpha_ch[chall[i]] + alpha_tg[targ[i]] +  X[i] * beta; 
  
  alpha ~ normal(0, 3);
  alpha_ch ~ normal(0, sigma_ch);
  sigma_ch ~ normal(0, 1); // half-normal
  alpha_tg ~ normal(0, sigma_tg);
  sigma_tg ~ normal(0, 1); // half-normal
  beta ~ student_t(5, 0, 2.5);
  
  y ~ bernoulli_logit(y_hat);
}
