library(rstan)
source("./models/models.r") 

# Find generated data then simulate from posterior distribution, if batch is true, batch data first
serial = function( model_name, dimension, data_dir ) {
    in_dir = paste0( data_dir, "datasets/", model_name, "/" )

    # Simulate from serial model
    x = read.table( paste0( in_dir, dimension ) )
    out_dir = paste0( data_dir, "samples/", model_name, "/", dimension, "/serial/" )
    dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )
    writeLines("Simulating from full model")
    hmc( x, model_name, out_dir, FALSE )
    return( 0 )
}

batch = function( model_name, dimension, seed, batch_num, data_dir ) {
    in_dir = paste0( data_dir, "datasets/", model_name, "/" )
    x = read.table( paste0( in_dir, dimension ) )
    set.seed(seed)
    
    # Randomly batch data
    N = nrow( x )
    n_batches = 20
    if ( N/n_batches - floor( N/n_batches ) > 10^(-4) ) {
        stop( "Wrong number of observations" )
    }
    batch_index = sample( rep( 1:n_batches, N/n_batches ), N )

    # Simulate from each batch of data
    out_dir = paste0( data_dir, "samples/",
                     model_name, "/", dimension, "/", seed, "/", batch_num, "/" )
    dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )
    batch_x = subset( x, batch_index == batch_num )
    writeLines( paste0( "Simulating from dimension ", dimension, 
                       " seed ", seed, " and batch ", batch_num ) )
    hmc( batch_x, model_name, out_dir, TRUE )
    return( 0 )
}

# Simulate from the posterior of dataset x under the specified model using hmc (STAN).
# Available models can be found in ./models/models.r
# batch determines if the sample has been batched and the prior needs to be flattened
hmc = function( x, model_name, out_dir, batch ) {

    # Setup model and output directories
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    n_iter = 10^4
    model = do.call( model_name, list( "x"=x, "batch"=batch ) )
    
    # Simulate from model using STAN
    mc_time = system.time( mc_out <- stan(file = paste0("./models/", model_name, ".stan"), 
              data = model, iter = 2*n_iter, chains = 1) )
    # Plot traceplot
    sample_plot = traceplot( mc_out )
    ggsave( filename = paste0( out_dir, "/traceplot.png" ), plot = sample_plot, dpi = 200 )

    # Save chain
    theta_samp = extract( mc_out )$theta
    if ( length( dim( theta_samp ) ) == 2 ) {
        write.table( theta_samp, file = paste0( out_dir, "theta" ), row.names=FALSE, col.names=FALSE )
    } 
    else if ( length( dim( theta_samp ) ) == 3 ) {
        write.table( theta_samp[,1,], file = paste0( out_dir, "theta-1" ), 
                    row.names = FALSE, col.names = FALSE )
        write.table( theta_samp[,2,], file = paste0( out_dir, "theta-2" ), 
                    row.names = FALSE, col.names = FALSE )
    }
    else {
        print( theta_samp )
        stop("Chain is the wrong number of dimensions, nothing saved")
    }
    write( mc_time[1], file = paste0( out_dir, "mc_time" ) )
    return( 0 )
}
