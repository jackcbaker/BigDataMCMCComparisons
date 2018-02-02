## Dependencies
library(mnormt)


## Constants
d = 2
out_dir = "/home/jbaker/documents/data/datasets/warped_gaussian/"
dir.create( out_dir, showWarnings=FALSE, recursive=TRUE )

theta = c( 0.5, 0 )
sigmay = 3
N = 800         # Number observations
start_seed = 13


## Simulate observations and store in DATADIR
set.seed( start_seed )
y = rnorm( N, theta[1] + theta[2]^2, sigmay )
write.table( y, file=paste0(out_dir, d), row.names=FALSE, col.names=FALSE )
