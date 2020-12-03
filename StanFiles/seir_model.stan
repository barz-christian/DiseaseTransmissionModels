//
// This Stan program defines a SEIR model
//
functions {
  real[] seir(real t, real[] y, real[] theta, 
              real[] x_r, int[] x_i) {
               
      // assign compartments and population size
      
      real S = y[1];
      real E = y[2];
      real I = y[3];
      real R = y[4];
      real N = x_i[1];
      
      // ODE parameters
      
      real beta = theta[1];
      real alpha = theta[2];
      real gamma = theta[3];
      
      // ODE system
      
      real dS_dt = -beta * I * S / N;
      real dE_dt =  beta * I * S / N - alpha * E;
      real dI_dt =  alpha * E - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days; // length observation period
  real y0[4];          // initial conditions of ODE
  real t0;             // inital time
  real ts[n_days];     // vector of time stamps
  int N;               // population size
  int cases[n_days];   // observation per time stamp
}

transformed data {
  real x_r[0];          // as there is nor real valued variable in the model it is left empty
  int x_i[1] = { N };   // the population size
}

parameters {
  real<lower=0> gamma;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(seir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
  //priors
  beta ~ normal(2, 1);
  alpha ~ uniform(0,21);
  gamma ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2), phi);
}
