include( "sg_methods.jl" )
include( "./models/models.jl" )

"""
Simulate from an LDA model using SGLD
"""
function sim_sgld( model::lda, subsize = 50::Int64, 
                   step_const = 1e-4::Float64, n_iter = 10^4::Int64 )

    # Generate objects and storage
    tuning = sgld( step_const, subsize )
    burn_in = round( Int, n_iter/10 )
    perplex_storage = zeros( ( floor( Int64, ( n_iter + burn_in )/100 ), 3 ) )

    # Simulate from model using SGLD
    tic()
    for ( tuning.iter in 1:( n_iter + burn_in ) )
        # Print progress, check the chain hasn't diverged, store perplexity
        if ( tuning.iter % 100 == 0 )
            perplex = perplexity( model, model.test )
            iter_time = toc()
            println("iter: $(tuning.iter) perplexity: $perplex time: $iter_time")
            perplex_storage[floor( Int64, tuning.iter/100 ),:] = [step_const,perplex,iter_time]
            tic()
            if ( sum( isnan( model.theta ) ) > 0 )
                print("\n")
                error("The chain has diverged")
            end
        end
        model = sgld_update( model, tuning )
    end

    return( perplex_storage )
end


"""
Simulate from an LDA model using SGHMC
"""
function sim_sghmc( model::lda, subsize = 50::Int64, 
                eta = 1e-4::Float64, alpha = 1.5::Float64, n_iter = 10^4::Int64 )

    # Generate objects and storage
    tuning = sghmc( subsize, alpha, eta )
    burn_in = round( Int, n_iter/10 )
    perplex_storage = zeros( ( floor( Int64, ( n_iter + burn_in )/100 ), 4 ) )

    # Simulate from model using SGLD
    tic()
    for ( tuning.iter in 1:( n_iter + burn_in ) )
        # Print progress, check the chain hasn't diverged, store perplexity
        if ( tuning.iter % 100 == 0 )
            perplex = perplexity( model, model.test )
            iter_time = toc()
            println("iter: $(tuning.iter) perplexity: $perplex time: $iter_time")
            open("./timings/time-$eta-$alpha", "a") do timefile
                write( timefile, "$(tuning.iter)\t$iter_time\t$perplex\t$eta\t$alpha\n" )
            end
            perplex_storage[floor( Int64, tuning.iter/100 ),:] = [eta,alpha,perplex,iter_time]
            tic()
            if ( sum( isnan( model.theta ) ) > 0 )
                print("\n")
                error("The chain has diverged")
            end
        end
        model = sghmc_update( model, tuning )
    end

    return( perplex_storage )
end

"""
Simulate from a model using SGLD
"""
function sim_sgld( model::statistical_model, 
                    step_const = 1e-4::Float64, n_iter = 10^4::Int64 )

    # Generate objects and storage
    tuning = sgld( step_const )
    burn_in = round( Int, n_iter/10 )
    output = zeros( ( ( n_iter + burn_in ), model.d ) )

    # Simulate from model using SGLD
    for ( tuning.iter in 1:( n_iter + burn_in ) )
        # Print progress, check the chain hasn't diverged
        if ( tuning.iter % 100 == 0 )
            print("$(tuning.iter) ")
            if ( sum( isnan( model.theta ) ) > 0 )
                print("\n")
                error("The chain has diverged")
            end
        end
        model = sgld_update( model, tuning )

        # Store new parameter value
        try
            output[tuning.iter,:] = model.theta
        catch LoadError
            output[tuning.iter,:] = model.theta[:,1]
        end
    end
    print("\n")

    # Discard burn-in
    output = slicedim( output, 1, burn_in:(n_iter + burn_in) )

    return( output )
end

"""
Simulate from a model using SGLD
"""
function sim_sghmc( model::statistical_model,
                eta = 1e-4::Float64, alpha = 2.0::Float64, n_iter = 10^4::Int64 )

    # Generate objects and storage
    tuning = sghmc( 200, alpha, eta )
    burn_in = round( Int, n_iter/10 )
    output = zeros( ( ( n_iter + burn_in ), model.d ) )

    # Simulate from model using SGHMC
    for ( tuning.iter in 1:( n_iter + burn_in ) )
        # Print progress
        if ( tuning.iter % 100 == 0 )
            print("$(tuning.iter) ")
            if ( sum( isnan( model.theta ) ) > 0 )
                print("\n")
                error("The chain has diverged")
            end
        end
        model = sghmc_update( model, tuning )

        # Store new parameter value
        try
            output[tuning.iter,:] = model.theta
        catch LoadError
            output[tuning.iter,:] = model.theta[:,1]
        end
    end
    print("\n")

    # Discard burn-in
    output = slicedim( output, 1, burn_in:(n_iter + burn_in) )

    return( output )
end
