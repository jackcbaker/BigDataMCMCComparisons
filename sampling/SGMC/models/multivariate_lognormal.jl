using Distributions

"""
Multivariate Lognormal with uninformative prior
"""
type multivariate_lognormal
    N::Int64                        # Number of observations
    d::Int64                        # Dimension
    x::Array{Float64,2}             # Data
    theta::Array{Float64,1}         # Location parameter
    sigma::Array{Float64,2}         # Scale parameter
    # Save computational time by storing regularly used transformations
    inv_sigma::Array{Float64,2}     # Inverse scale
end

"""
Simplified multivariate_lognormal declaration
"""
function multivariate_lognormal( x )
    (N, d) = size( x )
    theta = vec( rand( MvNormal( zeros(d), diagm( ones(d) ) ) ) )   # Simulate initial theta values
    sigma = diagm( ones( Float64, 2 ) )
    inv_sigma = inv( sigma )
    multivariate_lognormal( N, d, x, theta, sigma, inv_sigma )
end

"""
Gradient of log likelihood
"""
function dloglik( model::multivariate_lognormal, subsample )
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
function dlogdens( model::multivariate_lognormal, y )
    out = model.inv_sigma * ( log(y) - model.theta )
    return( out )
end

"""
Gradient of log prior
"""
function dlogprior( model::multivariate_lognormal )
    return( 0 )
end
