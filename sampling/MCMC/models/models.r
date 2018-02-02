# Build model files from data to pass to STAN. Use standard parameter values.

multivariate_gauss = function( x, batch ) {
    model = list( )
    model$N = nrow( x )
    model$d = ncol( x )
    model$x = x
    model$sigma = diag( rep( 1, model$d ) )
    return( model )
}

multivariate_t = function( x, batch ) {
    model = list( )
    model$N = nrow( x )
    model$d = ncol( x )
    model$x = x
    model$sigma = diag( c( 1, 1 ) )
    model$nu = 3
    return( model )
}

multivariate_lognormal = function( x, batch ) {
    model = list( )
    model$N = nrow( x )
    model$d = ncol( x )
    model$x = x
    model$sigma = diag( c( 1, 1 ) )
    return( model )
}

warped_gaussian = function( x, batch ) {
    model = list( )
    model$N = nrow( x )
    model$x = unlist( x )
    names( model$x ) = NULL
    model$sigmay = 3
    model$prior_mean = c( 0.0, 0 )
    model$prior_scale = diag( c( 0.1, 0.1 ) )
    # Flatten prior if it's a subposterior
    if ( batch == TRUE ) {
        model$prior_scale = model$prior_scale*20
    }
    return( model )
}

correlated_mixture = function( x, batch ) {
    model = list( )
    model$K = 2
    model$d = ncol( x )
    model$N = nrow( x )
    model$x = x
    model$sigma = array( rep( NA, model$K*model$d*model$d ), c(model$K,model$d,model$d) )
    model$sigma[1,,] = matrix( c( 5,3,3,5 ), ncol=2 )
    model$sigma[2,,] = matrix( c( 5,-3,-3,5 ), ncol=2 )
    model$theta0 = rep( 0, model$d )
    model$theta_scale0 = diag( rep( 100, model$d ) )
    # Flatten prior if it's a subposterior
    if ( batch == TRUE ) {
        model$mu_scale0 = model$mu_scale0*20
    }
    return( model )
}

gauss_mixture = function( x, batch ) {
    model = list( )
    model$K = 2
    model$d = ncol( x )
    model$N = nrow( x )
    model$x = x
    model$sigma = array( rep( NA, model$K*model$d*model$d ), c(model$K,model$d,model$d) )
    model$sigma[1,,] = diag( rep( 1, model$d ) )
    model$sigma[2,,] = diag( rep( 1, model$d ) )
    model$theta0 = rep( 0, model$d )
    model$theta_scale0 = diag( rep( 100, model$d ) )
    # Flatten prior if it's a subposterior
    if ( batch == TRUE ) {
        model$mu_scale0 = model$mu_scale0*20
    }
    return( model )
}
