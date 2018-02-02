data {
    int<lower=1> N;         // Number of observations
    int<lower=1> d;         // Number of dimensions
    row_vector[d] x[N];     // Observations
    matrix[d,d] sigma;      // Scale
}
parameters {
    row_vector[d] theta;    // Location
}
model {
    x ~ multi_normal( theta, sigma );
}
