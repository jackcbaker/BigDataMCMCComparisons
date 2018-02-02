using Distributions

"""
Multivariate t with uninformative prior
"""
type multivariate_t
    N::Int64                        # Number observations
    d::Int64                        # Dimension
    x::Array{Float64,2}             # Data
    theta::Array{Float64,1}         # Location parameter
    sigma::Array{Float64,2}         # Scale parameter
    nu::Float64                     # Degrees of freedom
    # Save computational time by storing regularly used transformations
    inv_sigma::Array{Float64,2}     # Inverse scale
end

"""
Simplified multivariate_t declaration
"""
function multivariate_t( x )
    N = size( x, 1 )
    d = size( x, 2 )
    theta = vec( rand( MvNormal( zeros(d), 10*eye(d) ) ) )  # Simulate initial theta vals
    sigma = diagm( ones( Float64, 2 ) )
    nu = 3
    inv_sigma = inv( sigma )
    multivariate_t( N, d, x, theta, sigma, nu, inv_sigma )
end

"""
Gradient of log likelihood
"""
function dloglik( model::multivariate_t, subsample )
    out = zeros( Float64, model.d )
    for ( i in subsample )
        y = vec( model.x[i,:] )
        out += vec( dlogdens( model, y ) )
    end
    return( out )
end

"""
Gradient of log density
"""
function dlogdens( model::multivariate_t, y )
    out = ( model.nu + model.d ) / ( model.nu ) * model.inv_sigma * ( y - model.theta )
    out /= 1 + 1/model.nu * ( y - model.theta )' * model.inv_sigma * ( y - model.theta )
    return( out )
end

"""
Gradient of log prior
"""
function dlogprior( model::multivariate_t )
    return( 0 )
end
