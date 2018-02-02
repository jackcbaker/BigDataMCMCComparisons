source("parallel_methods.r")
source("diagnostics.r")

tuning = function( model_name, dimension ) {

    # Load in data
    seed = 2
    in_dir = paste0( "/home/jbaker/documents/data/comparisons/samples/",
                     model_name, "/", dimension, "/", seed, "/" )
    if ( model_name == "gaussian_mixture" ) {
        recomb_mixture( dimension, seed, method_name )
    }
    batch_array = array( rep( NA, 10^4*dimension*20 ), c( 10^4, dimension, 20 ) )
    for ( batch_num in 1:20 ) {
        batch_array[,,batch_num] = as.matrix( read.table( paste0( in_dir, batch_num, "/theta" ) ) )
    }

    # Recombine using kdemc
    out_dir = paste0( "/home/jbaker/documents/data/comparisons/tuning/",
                     model_name, "/", dimension, "/kdemc/" )
    dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )
    for ( bw in c( 0.1, 0.2, 0.3, 0.4, 0.5 ) ) {
        writeLines( paste0( "Combining with bandwidth ", bw ) )
        sample = kdemc( batch_array, bw*diag( rep( 1, dimension ) ) )
        write.table( sample, paste0( out_dir, bw ), row.names=FALSE, col.names=FALSE )
    }
    tuning_diagnostics( "multivariate_t", "kdemc" )
    return(0)
}

recombine = function( model_name, dimension, bw ) {

    for ( seed in 1:10 ) {
        # Load in data
        in_dir = paste0( "/home/jbaker/documents/data/comparisons/samples/",
                         model_name, "/", dimension, "/", seed, "/" )
        if ( model_name == "gaussian_mixture" ) {
            recomb_mixture( dimension, seed, method_name )
        }
        batch_array = array( rep( NA, 10^4*dimension*20 ), c( 10^4, dimension, 20 ) )
        for ( batch_num in 1:20 ) {
            batch_array[,,batch_num] = as.matrix( read.table( paste0( in_dir, batch_num, "/theta" ) ) )
        }

        # Combine using Consensus MC
        out_dir = paste0( "/home/jbaker/documents/data/comparisons/samples/",
                         model_name, "/", dimension, "/", "consensus/", seed, "/" )
        dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )
        sample = consensus( batch_array )
        write.table( sample, paste0( out_dir, "theta" ), row.names=FALSE, col.names=FALSE )
        # Combine using KDEMC
        out_dir = paste0( "/home/jbaker/documents/data/comparisons/samples/",
                         model_name, "/", dimension, "/", "kdemc/", seed, "/" )
        dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )
        sample = kdemc( batch_array, bw*diag( rep( 1, dimension ) ) )
        write.table( sample, paste0( out_dir, "theta" ), row.names=FALSE, col.names=FALSE )
    }
}
