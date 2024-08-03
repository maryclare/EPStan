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
  vector[n2] z2; // Regression coefficients
  vector[n2] lxi; // log scales
  vector[n2] logitdelta; // logit zolotorev random variables 
}
transformed parameters {
  vector[n2] xi = exp(lxi);
  vector[n2] delta = pi()./(1.0 + exp(-logitdelta));
}
model {
  y ~ normal(x*z2, sigma); 
  target += -n2*log(2.0*pi());
  target += n2*(log(lambda)/q + log(q) - lgamma(1.0/q)); 
  target += -sum(xi);
  target += -lambda^(2.0/q)*sum(z2^2.0 .* sin(q*delta/2.0) .* sin((2.0 - q)*delta/2.0)^((2.0 - q)/q) .* sin(delta)^(-2.0/q) .* xi^((q - 2.0)/q));
  target += sum(lxi) + sum(logitdelta) - 2.0*sum(log((1.0 + exp(logitdelta)))); // Change of variables for conversion to reals
}
