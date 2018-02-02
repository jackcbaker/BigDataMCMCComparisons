library(mnormt)

# Simulate data from a multivariate t-distribution, with the desired parameter estimates
set.seed(13)
out_dir = "/home/jbaker/documents/data/comparisons/datasets/multivariate_t/"
dir.create( out_dir, showWarnings = FALSE, recursive = TRUE )

dataset = rmt( 800, mean = c( 0, 0 ), S = diag( c( 1, 1 ) ), df = 3 )
write.table( dataset, paste0( out_dir, 2 ), row.names=FALSE, col.names=FALSE )
