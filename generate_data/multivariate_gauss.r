library(mnormt)

# Simulate data from a multivariate t-distribution, with the desired parameter estimates
set.seed(13)
out_dir = "/home/jbaker/documents/data/comparisons/datasets/multivariate_gauss/"
dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )

for ( dimension in 2:10 ) {
    dataset = rmnorm( 800, mean = rep( 0, dimension ), varcov = diag( rep( 1, dimension ) ) )
    write.table( dataset, paste0( out_dir, dimension ), row.names=FALSE, col.names=FALSE )
}
