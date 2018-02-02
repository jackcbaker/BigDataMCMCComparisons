using Distributions

"""
Gaussian mixture model with common scale and Guassian prior
"""
type gauss_mixture
    N::Int64                        # Number of observations
    d::Int64                        # Dimension
    K::Int64                        # Number of mixture components
    x::Array{Float64,2}             # Data
    theta::Array{Float64,2}         # Location parameters, columns are components, rows dimension
    sigma::Array{Float64,2}         # Scale parameter (common)
    theta0::Array{Float64,1}        # Location hyperparameter
    sigma0::Array{Float64,2}        # Scale hyperparameter
    # Save computational time by storing regularly used transformations
    inv_sigma::Array{Float64,2}     # Inverse scale
    inv_sigma0::Array{Float64,2}    # Inverse prior scale
end

"""
Set some default values
"""
function gauss_mixture( x )
    ( N, d ) = size( x )
    K = 2
    sigma = diagm( ones( Float64, 2 ) )
    sigma0 = 100*eye(d)
    theta0 = vec( zeros( Float64, 2 ) )
    theta = rand( MvNormal( theta0, eye(d) ), K )    # Simulate initial values from a guess
    theta = [0.1 0; 0.1 0]
    gauss_mixture( N, d, K, x, theta, sigma, theta0, sigma0 )
end

"""
Calculate transformations
"""
function gauss_mixture( N, d, K, x, theta, sigma, theta0, sigma0 )
    inv_sigma = inv( sigma )
    inv_sigma0 = inv( sigma0 )
    gauss_mixture( N, d, K, x, theta, sigma, theta0, sigma0, inv_sigma, inv_sigma0 )
end

"""
Gradient of log likelihood
"""
function dloglik( model::gauss_mixture, subsample )
    out = zeros( Float64, ( model.d, model.K ) )
    for ( i in subsample )
        y = model.x[i,:]
        out += dlogdens( model, y )
    end
    return( out )
end

"""
Gradient of log density
"""
function dlogdens( model::gauss_mixture, y )
    out = zeros( Float64, ( model.d, model.K ) )
    denom = 0
    for ( k in 1:model.K )
        out[:,k] = model.inv_sigma * ( vec( y ) - vec( model.theta[:,k] ) ) * 
            pdf( MvNormal( vec( model.theta[:,k] ), model.sigma ), vec( y ) )
        denom += pdf( MvNormal( vec( model.theta[:,k] ), model.sigma ), vec( y ) )
    end
    out /= denom
    return( out )
end

"""
Gradient of log prior
"""
function dlogprior( model::gauss_mixture )
    out = zeros( Float64, ( model.d, model.K ) )
    for ( k in 1:model.K )
        out[:,k] = model.inv_sigma0 * ( vec( model.theta[:,k] ) - model.theta0 )
    end
    return( out )
end

"""
Score variance estimate
"""
function grad_ests( model::gauss_mixture, subsample )
    score_var = zeros( ( model.d, model.d, model.K ) )
    subsize = length( subsample )
    grad_matrix = gradient_matrix( model, subsample )
    grad_ave = sum( gradient_matrix, 1 )/subsize
    grad_matrix = grad_matrix - grad_ave
    for ( i in 1:subsize )
        for ( k in 1:model.K )
            score_var[:,:,k] += vec( grad_matrix[i,:,k] ) * vec( grad_matrix[i,:,k] )'
        end
    end
    score_var /= subsize - 1
    return( score_var )
end

"""
Calculate Gradient Matrix
"""
function gradient_matrix( model::gauss_mixture, subsample )
    subsize = length( subsample )
    out = zeros( Float64, (subsize, model.d, model.K ) )
    for ( i in subsample )
        y = model.x[i,:]
        out[i,:,:] = dlogdens( model, y )
    end
    return( out )
end
