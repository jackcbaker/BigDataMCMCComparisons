include("multivariate_gauss.jl")
include("multivariate_t.jl")
include("multivariate_lognormal.jl")
include("gauss_mixture.jl")
include("warped_gaussian.jl")
include("lda.jl")

statistical_model = Union{ multivariate_gauss, multivariate_t,  gauss_mixture, multivariate_lognormal, warped_gaussian }
