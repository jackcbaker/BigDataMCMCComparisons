data {
    int<lower=1> K;         // number of mixture components
    int<lower=1> d;         // dimension
    int<lower=1> N;         // number of data points
    row_vector[d] x[N];     // observations
    matrix[d,d] sigma[K];   // scales of mixture components (common)
    vector[d] theta0;          // hyperparameters for location
    matrix[d,d] theta_scale0;  // Scale hyperparameter for location
}
parameters {
    row_vector[d] theta[K];    // locations of mixture components
}
model {
    real ps[K];             // temp for log component densities
    for (k in 1:K) {
        theta[k] ~ multi_normal( theta0, theta_scale0 );
    }
    for (n in 1:N) {
        for (k in 1:K) {
            ps[k] <- log(1.0/K) + multi_normal_log(x[n],theta[k],sigma[k]);
        }
    increment_log_prob(log_sum_exp(ps));
    }
}
