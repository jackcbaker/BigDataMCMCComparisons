## Dependencies
library(mnormt)


## Constants
K = 2           # Num mixture components
d = 2           # Dimension
mu = matrix( rep(NA, d*K), ncol=d )   # Matrix of means
out_dir = "/home/jbaker/documents/data/comparisons/datasets/gauss_mixture/"
dir.create( out_dir, showWarnings=FALSE, recursive=TRUE )

mu[1,] = c(0.1,0.1)
mu[2,] = c(0,0)
sigma = diag(rep(1,d))       # Common variance
N = 10000         # Number observations
s = c(0.5,0.5)  # Label
start_seed = 13


## Simulate labels, observations and store in DATADIR
set.seed( start_seed )
z = sample( 1:K, N, prob=s, replace=TRUE )
y = t( apply( mu[z,], 1, function (x) rmnorm( 1, x, sigma ) ) )
write.table( y, file=paste0(out_dir, d), row.names=FALSE, col.names=FALSE )
