# Dependencies
library(FNN)

# Obtain the KL divergence from the true posterior for the given method
kl_divergence = function( model, dimension, method, data_dir ) {

    # Read in 'exact' data
    writeLines(method)
    true_sample = paste0( data_dir, "samples/", model, "/", dimension, "/serial/theta" )
    true_sample = as.matrix( read.table( true_sample ) )

    # Read in approx data, calculate KL divergence from exact data, store
    approx_sample = matrix( rep( NA, dimension*10^4 ), 10^4, dimension )
    KL_estimates = rep( NA, 10 )
    approx_dir = paste0( data_dir, "samples/", model, "/", dimension, "/", method, "/" )
    for ( batch in 1:10 ) {
        approx_sample = as.matrix( read.table( paste0( approx_dir, batch, "/theta" ) ) )
        KL_estimates[batch] = KL.divergence( true_sample, approx_sample, k=100 )[100]
        if ( KL_estimates[batch] < 0 ) {
            writeLines( "KL Divergence estimate is negative" )
            KL_estimates[batch] = NA
        }
    }
    print( KL_estimates )

    return( KL_estimates )
}

# Plot KL divergence plots for each method
plot_kl = function( model, dimension, data_dir ) {

    # Obtain KL estimates for each method
    methods = c( "consensus", "kdemc", "sgld", "sghmc" )
    kl_frame = data.frame( "method" = rep( methods, each=10 ), 
                           "kl_divergence" = rep( NA, 10*length(methods) ) )
    for ( method in methods ) {
        kl_estimates = kl_divergence( model, dimension, method, data_dir )
        kl_frame$kl_divergence[ kl_frame$method == method ] = kl_estimates
    }

    # Create a boxplot from the results
    out_path = paste0( data_dir, "plots/", model, "/", dimension, "/" )
    dir.create( out_path, recursive = TRUE, showWarnings = FALSE )
    pdf( paste0( out_path, "boxplot.pdf" ), width=8, height=4 )
        boxplot( kl_frame$kl_divergence~kl_frame$method, 
                 col = c( "yellow", "yellow", "green", "green" ),
                 xlab = "Method", ylab = "KL-divergence" )
        legend( "topright", c( "Parallel", "Stochastic Gradient" ), col=c("yellow", "green"),
                pch=16 )
    dev.off()
}
