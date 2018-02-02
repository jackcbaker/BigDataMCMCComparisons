# Dependencies
library( ggplot2 )

# Stack data frames for each method on top of one another for the final plot
stack_frames = function( model, dimension, data_dir ) {
    in_dir = paste0( data_dir, "samples/", model, "/", dimension, "/" )
    methods = c( "consensus", "kdemc", "sgld", "sghmc" )
    serial_frame = read.table( paste0( in_dir, "serial/theta" ) )
    colnames( serial_frame ) = paste0( "theta", 1:2 )
    serial_frame$method = "serial"
    serial_frame$type = "Exact"
    
    full_frame = data.frame( "theta1" = {}, "theta2" = {}, "method" = {}, "type"={} )
    for ( method in methods ) {
        current_frame = read.table( paste0( in_dir, method, "/2/theta" ) )
        colnames( current_frame ) = paste0( "theta", 1:2 )
        current_frame$method = method
        if ( method == "consensus" || method == "kdemc" ) {
            current_frame$type = "Parallel"
        }
        else {
            current_frame$type = "Stochastic Gradient"
        }
        serial_frame$method = method
        full_frame = rbind( full_frame, current_frame, serial_frame )
    }

    return( full_frame )
}

# Plot empirical contours for each method
contour_comparison = function( stacked_frame, model, dimension, data_dir ) {
    comparison_plot = ggplot( stacked_frame, aes( x=theta1, y=theta2, color=type ) ) +
        stat_density2d( alpha=0.6, size=1.5 ) +
        facet_grid( . ~ method ) +
        scale_colour_manual( values = c( "red", "blue", "green" ) )
    dir.create(paste0( data_dir, "plots/",model,"/",dimension,"/"),
               recursive = TRUE, showWarnings = FALSE )
    ggsave(paste0( data_dir, "plots/", model, "/", dimension, "/", "contours.pdf"),
           comparison_plot, width = 8, height = 2.5)
    return( comparison_plot )
}

# Contour plots
comparison_plot = function( model, dimension, data_dir ) {
    plot_frame = stack_frames( model, dimension, data_dir )
    out_plot = contour_comparison( plot_frame, model, dimension, data_dir )
    return( out_plot )
}
