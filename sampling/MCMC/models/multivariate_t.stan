data {
    int<lower=1> N;         // Number of observations
    int<lower=1> d;         // Number of dimensions
    row_vector[d] x[N];     // Observations
    matrix[d,d] sigma;      // Scale
    real<lower=0> nu;       // DF
}
parameters {
    row_vector[d] theta;    // Location
}
model {
    x ~ multi_student_t( nu, theta, sigma );
}
