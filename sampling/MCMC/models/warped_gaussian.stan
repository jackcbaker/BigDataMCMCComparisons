data {
    int<lower=1> N;             // Number of Observations
    real x[N];            // Observations
    real<lower=0> sigmay;       // Scale of y
    row_vector[2] prior_mean;            // Mean of prior for theta
    matrix[2,2] prior_scale;     // Scale of theta
}
parameters {
    row_vector[2] theta;        //location
}
model {
    theta ~ multi_normal( prior_mean, prior_scale );
    x ~ normal( theta[1] + theta[2]^2, sigmay );
}
