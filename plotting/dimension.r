# Dependencies
library(FNN)
library(ggplot2)

# Obtain the KL divergence from the true posterior for the given method
kl_divergence = function( model, dimension, method, data_dir ) {

    # Read in 'exact' data
    writeLines( paste0( "method: ", method, ", dimension: ", dimension ) )
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

    return( KL_estimates )
}

# Plot KL divergence plots for each method
plot_kl = function( model, data_dir ) {

    # Obtain KL estimates for each method
    methods = c( "consensus", "kdemc", "sgld", "sghmc" )
    dimensions = 2:10
    seeds = 1:10
    kl_frame = expand.grid( methods, dimensions )
    colnames( kl_frame ) = c( "method", "dimension" )
    kl_frame$kl = NA
    kl_frame$kl_lower = NA
    kl_frame$kl_upper = NA
    for ( current_method in methods ) {
        for ( current_dimension in dimensions ) {
            kl_estimates = kl_divergence( model, current_dimension, current_method, data_dir )
            index = with( kl_frame, method == current_method & dimension == current_dimension )
            kl_frame$kl[index] = mean( kl_estimates, na.rm=TRUE )
            kl_frame$kl_lower[index] = quantile( kl_estimates, 0.05, na.rm = TRUE )
            kl_frame$kl_upper[index] = quantile( kl_estimates, 0.95, na.rm = TRUE )
        }
    }

    # Create a ribbon plot from the results
    out_dir = paste0( data_dir, "plots/", model, "/" )
    dir.create( out_dir, recursive = TRUE, showWarnings = FALSE )
    out_path = paste0( out_dir, "dimension.pdf" )
    ggplot( kl_frame, aes( x = dimension, y = kl, fill = method, color = method ) ) +
        geom_line() +
        geom_ribbon( aes( ymin = kl_lower, ymax = kl_upper ), alpha = 0.2 ) +
        xlab("Dimension") +
        ylab("KL-divergence")
    ggsave( out_path, width = 8, height = 4 )
}
