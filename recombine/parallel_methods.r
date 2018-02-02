library(parallelMCMCcombine)
library(mnormt)

consensus = function( batch_container ) {

    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )
    s = dim( batch_container )[3]

    temp_container = array( rep( NA, s*p0*sample_size ), c(p0,sample_size,s) )
    for ( i in 1:s ) {
        temp_container[,,i] = t( batch_container[,,i] )
    }

    sample_consens = consensusMCindep( temp_container )

    return( t( sample_consens ) )
}

kdemc = function( batch_container, bw_init ) {

    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    s = dim( batch_container )[3]
    sample_size = nrow( batch_container[,,1] )
    # Ensure bw_init is diagonal
    bw_init = diag( diag( bw_init ) )
    
    ## Initialise
    sample_container = matrix( rep( NA, p0*sample_size ), ncol=p0 )
    T_current = sample( sample_size, s, replace=TRUE )
    bw = bw_init
    w_t = calc_weights( T_current, batch_container, s, bw, p0 )
    for ( j in 1:sample_size ) {
        print_progress( j )
        bw = ( j^(-1/(4+p0)) )*bw_init
        for ( i in 1:s ) {
            C_current = T_current
            C_current[i] = sample( sample_size, 1 )
            u = runif(1)
            w_c = calc_weights( C_current, batch_container, s, bw, p0  )
            if ( log(u) < w_c - w_t ) {
                T_current = C_current
                w_t = w_c
            }
        }
        mean_current = calc_mean( T_current, batch_container, s, bw )
        sample_container[j,] = rmnorm( 1, mean=mean_current, var=bw^2/s )
    }
    cat("\n")

    return( sample_container )
}

# Helper functions for recomb_npara
calc_weights = function( t_lab, batch_container, s, bw, p0 ) {
    theta_vals = matrix( rep( NA, s*p0 ), ncol=p0 )
    for ( i in 1:s ) {
        theta_vals[i,] = batch_container[t_lab[i],,i]
    }
    theta_mean = colMeans( theta_vals )
    # Store log weights for numerical stability
    weight_val = sum( dmnorm( theta_vals, mean = theta_mean, var = bw^2, log = TRUE ) )

    return( weight_val )
}

calc_mean = function( t_lab, batch_container, s, bw ) {
    theta_mean = 0
    for ( i in 1:s ) {
        theta_mean = theta_mean + batch_container[t_lab[i],,i]/s
    }

    return( theta_mean )
}

print_progress = function( j ) {
    if ( j %% 1000 == 0 ) {
        cat( paste0( j, " " ) )
    }
}
