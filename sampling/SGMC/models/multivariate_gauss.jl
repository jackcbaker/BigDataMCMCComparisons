using Distributions

"""
Multivariate Gaussian with uninformative prior
"""
type multivariate_gauss
    N::Int64                        # Number of observations
    d::Int64                        # Dimension
    x::Array{Float64,2}             # Data
    theta::Array{Float64,1}         # Location parameter
    sigma::Array{Float64,2}         # Scale parameter
    # Save computational time by storing regularly used transformations
    inv_sigma::Array{Float64,2}     # Inverse scale
end

"""
Simplified multivariate_gauss declaration
"""
function multivariate_gauss( x )
    (N, d) = size( x )
    theta = rand( MvNormal( zeros(d), 10*eye(d) ) )   # Simulate initial theta values
    sigma = eye(d)
    inv_sigma = inv( sigma )
    multivariate_gauss( N, d, x, theta, sigma, inv_sigma )
end

"""
Gradient of log likelihood
"""
function dloglik( model::multivariate_gauss, subsample )
    out = zeros( Float64, model.d )
    for ( i in subsample )
        y = vec( model.x[i,:] )
        out += dlogdens( model, y )
    end
    return( out )
end

"""
Gradient of log density
"""
function dlogdens( model::multivariate_gauss, y )
    out = model.inv_sigma * ( y - model.theta )
    return( out )
end

"""
Gradient of log prior
"""
function dlogprior( model::multivariate_gauss )
    return( 0 )
end
