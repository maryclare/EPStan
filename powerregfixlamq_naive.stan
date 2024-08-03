data {
  int<lower=0> m;   // Number of observations
  int<lower=0> n2;   // Number of regression coefficients
  vector[m] y;       // Response
  matrix[m, n2] x;    // Predictor(s)
  real<lower=0> sigma;
  real<lower=0> lambda; // rate
  real<lower=0, upper = 2.0> q; //Shape parameter
}
parameters {
  vector[n2] w; // Regression coefficients
}
transformed parameters {
  vector[n2] z2 = w * lambda^(-1.0/q);
}
model {
  y ~ normal(x*z2, sigma); 
  for (j in 1:n2) {
    target += log(q) - abs(w[j])^q - log(2.0) - lgamma(1.0/q);
  }
}
