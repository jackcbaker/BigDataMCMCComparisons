using Distributions
include("./models/models.jl")

"""
Container for sgld parameters
"""
type sgld
    subsize::Int64                  # Minibatch size
    step_const::Float64             # Stepsize tuning constant
    iter::Int64
    subsample::Array{ Int64, 1 }
end

"""
Set standard tuning values
"""
function sgld( step_const, subsize = 200::Int64 )
    iter = 1
    subsample = [-1]
    sgld( subsize, step_const, iter, subsample )
end

"""
Container for sghmc parameters
"""
type sghmc
    subsize::Int64           # Minibatch size
    alpha::Float64                  # Momentum param
    eta::Float64                    # Learning param
    L::Int64                        # Trajectory
    iter::Int64
    subsample::Array{ Int64, 1 }
end

"""
Set standard tuning values for sghmc
"""
function sghmc( subsize, alpha, eta )
    iter = 1
    subsample = [-1]
    L = 3
    sghmc( subsize, alpha, eta, L, iter, subsample )
end

"""
Calculate stepsize for sgld
"""
function stepsize( tuning::sgld )
    tuning.step_const*( 1 + tuning.iter )^(-0.33)
end

"""
Update one step of Stochastic Gradient Langevin Dynamics
"""
function sgld_update( model::statistical_model, tuning::sgld )
    
    # Subsample data
    tuning.subsample = sample( 1:model.N, tuning.subsize )
    println( tuning.subsize )
    error()

    # Simulate Langevin dynamics of new submodel
#    epsilon = stepsize( tuning )
    epsilon = tuning.step_const
    injected_noise = MvNormal( zeros(model.d), eye(model.d) )
    # Update parameter
    model.theta += epsilon/2 * dlogprior( model )
    model.theta += epsilon/2 * ( model.N/tuning.subsize ) * dloglik( model, tuning.subsample )
    try
        model.theta += sqrt( epsilon ) * rand( injected_noise )
    catch LoadError
        model.theta += rand( injected_noise, size( model.theta, 2 ) )
    end

    return( model )
end

"""
Update one step of Stochastic Gradient Hamiltonian Monte Carlo
"""
function sghmc_update( model::statistical_model, tuning::sghmc )
    
    # Subsample data
    tuning.subsample = sample( 1:model.N, tuning.subsize )
    println( tuning.subsize )
    error()

    # Generate new momentum values and reparameterize
#    if ( size( model.theta, 2 ) == 1 )
#        v = rand( MvNormal(eye(model.d)) )
#    else
#        v = rand( MvNormal(eye(model.d)), size( model.theta, 2 ) )
#    end
#    v *= sqrt( tuning.eta )
#    v = zeros( size( model.theta ) )
    v = sqrt( tuning.eta ) * rand( Normal(), size( model.theta ) )
#    noise_param = tuning.alpha*tuning.eta/4

    # Simulate from SGHMC dynamics
    for ( i = 1:tuning.L )
        model.theta += v
        # Calculate potential energy
        PE = model.N/tuning.subsize * dloglik( model, tuning.subsample ) + dlogprior( model )
        # Momentum step
#        noise_distn = MvNormal( zeros(model.d), noise_param * eye(model.d) )
        injected_noise = sqrt( 2 * tuning.alpha * tuning.eta ) * rand( Normal(), size(model.theta) )
        v += tuning.eta*PE - tuning.alpha*v + injected_noise
    end

    return( model )
end

"""
Update one step of Stochastic Gradient Langevin Dynamics for a Latent Dirichlet Allocation model
"""
function sgld_update( model::lda, tuning::sgld )
    
    # Subsample documents
    tuning.subsample = sample( 1:model.M, tuning.subsize )

    # Simulate Langevin dynamics of new submodel
    epsilon = stepsize( tuning )
    injected_noise = Normal( 0, epsilon )
    model.zcounts[tuning.subsample,:,:] = update_topics( model, tuning.subsample )
    # Topic word parameter
    model.theta += epsilon/2 * grad_theta_prior( model )
    model.theta += epsilon/2 * ( model.M/tuning.subsize ) * grad_theta( model, tuning.subsample )
    model.theta += rand( injected_noise, ( model.K, model.V ) )

    return( model )
end

"""
Update one step of Stochastic Gradient Hamiltonian MC for a Latent Dirichlet Allocation model
"""
function sghmc_update( model::lda, tuning::sghmc )
    
    # Subsample data
    tuning.subsample = sample( 1:model.M, tuning.subsize )

    # Generate new momentum values and reparameterize
    v = rand( Normal( ), size( model.theta ) )
    v *= sqrt( tuning.eta )
    noise_param = tuning.alpha*tuning.eta/4
    model.zcounts[tuning.subsample,:,:] = update_topics( model, tuning.subsample )

    # Simulate from SGHMC dynamics
    for ( i = 1:tuning.L )
        model.theta += v
        # Calculate potential energy
        PE = grad_theta( model, tuning.subsample )
        PE *= model.M/tuning.subsize
        PE += grad_theta_prior( model )
        # Momentum step
        v += tuning.eta*PE - tuning.alpha*v + rand( Normal( 0, noise_param ), size( model.theta ) )
    end

    return( model )
end
