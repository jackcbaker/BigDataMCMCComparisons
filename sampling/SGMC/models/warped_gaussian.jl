using Distributions

"""
Multivariate Gaussian with uninformative prior
"""
type warped_gaussian
    N::Int64                        # Number of observations
    d::Int64                        # Dimension
    x::Array{Float64,2}             # Data
    theta::Array{Float64,1}         # Location parameter
    sigmay::Float64                 # Scale of y
    prior_mean::Array{Float64,1}    # Prior mean for theta
    prior_scale::Array{Float64,2}            # Prior scale for theta
    inv_prior_scale::Array{Float64,2}        # Store inverse prior scale for efficiency     
end

"""
Simplified warped_gaussian declaration
"""
function warped_gaussian( x )
    N = size( x, 1 )
    d = 2
    prior_mean = [0.0,0.0]
    prior_scale = diagm( [0.1,0.1] )
    inv_prior_scale = inv( prior_scale )
    sigmay = 3
    theta = rand( MvNormal( [0.5,0.0], eye(d) ) )   # Simulate initial theta values
    warped_gaussian( N, d, x, theta, sigmay, prior_mean, prior_scale, inv_prior_scale )
end

"""
Gradient of log likelihood
"""
function dloglik( model::warped_gaussian, subsample )
    out = zeros( Float64, model.d )
    for ( i in subsample )
        y = model.x[i]
        out += dlogdens( model, y )
    end
    return( out )
end

"""
Gradient of log density
"""
function dlogdens( model::warped_gaussian, y )
    out = zeros( model.d )
    out[1] = 1/model.sigmay^2 * ( y - model.theta[1] - model.theta[2]^2 )
    out[2] = 2*model.theta[2]/model.sigmay^2 * ( y - model.theta[1] - model.theta[2]^2 )
    return( out )
end

"""
Gradient of log prior
"""
function dlogprior( model::warped_gaussian )
    out = zeros( model.d )
    out = - model.inv_prior_scale * ( model.theta - model.prior_mean )
    return( out )
end

"""
Potential energy
"""
function pe( model::warped_gaussian, current_theta )
    out = 0
    for ( i in 1:model.N )
        y = model.x[i]
        out += logdens( model, current_theta, y )
    end
    out += logprior( model, current_theta )

    return( -out )
end

function logdens( model::warped_gaussian, current_theta, y )
    out = - 1/(2*model.sigmay^2) * ( y - current_theta[1] - current_theta[2]^2 )^2
    return( out )
end

function logprior( model::warped_gaussian, current_theta )
    out = ( current_theta - model.prior_mean )' * 
            model.inv_prior_scale * ( current_theta - model.prior_mean )
    return( out )
end
