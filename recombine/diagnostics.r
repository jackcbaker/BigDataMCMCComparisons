library(FNN)
library(ggplot2)

tuning_diagnostics = function( model_name, method ) {
    cat("Creating tuning diagnostics\nparameter: ")
    in_dir = paste0("/home/jbaker/documents/data/comparisons/tuning/", model_name, "/2/", method)
    stacked = read.table( paste0( "/home/jbaker/documents/data/comparisons/samples/multivariate_t/2/serial/theta" ) )
    truth = as.matrix( stacked )
    colnames(stacked) = paste0( "theta", 1:2 )
    stacked$bw = 0
    params = dir( in_dir, pattern = "[[:digit:]]" )
    kl_frame = data.frame( params = params, kl = rep( NA, length( params ) ) )
    for ( param in params ) {
        cat(paste0(param," "))
        sample_curr = read.table( paste0( in_dir, "/", param ) )
        kl_estimate = KL.divergence( truth, as.matrix( sample_curr ), k=10^2 )
        kl_estimate = kl_estimate[100]
        kl_frame[ params == param, 2 ] = kl_estimate
        colnames(sample_curr) = paste0( "theta", 1:2 )
        sample_curr$bw = param
        stacked = rbind( stacked, sample_curr )
    }
    p = ggplot( stacked, aes( x=theta1, y=theta2 ) ) +
        geom_density2d() +
        facet_grid( bw~. )
    ggsave( paste0( in_dir, "/tune_plot.pdf" ) )
    write.table( kl_frame, paste0( in_dir, "/kl.dat" ) )
}
